#include "main.h"

int writeCigar(char** o_buf, int* o_buflen, int count, char code,enum CigarFormat format)
{
    int i = 0;
    if (count <= 0)
    {
        return 1;
    }
    if (format == EXPANDED_CIGAR_STRING)
    {
        int n = min(*o_buflen, count);
        for (i = 0; i < n; i++)
        {
            *(*o_buf)++ = code;
        }
        *o_buflen -= n;
        if (*o_buflen == 0)
        {
            *(*o_buf - 1) = '\0';
        }
        return *o_buflen > 0;
    }
    else if (format == COMPACT_CIGAR_STRING)
    {
        if (*o_buflen == 0)
        {
            *(*o_buf - 1) = '\0';
            return 0;
        }
        int written = snprintf(*o_buf, *o_buflen, "%d%c", count, code);
        if (written > *o_buflen - 1)
        {
            *o_buf = '\0';
            return 0;
        }
        else
        {
            *o_buf += written;
            *o_buflen -= written;
            return 1;
        }
    }
    else if (format == COMPACT_CIGAR_BINARY)
    {
        // binary format with non-zero count byte followed by char (easier to examine programmatically)
        while (1)
        {
            if (*o_buflen < 3)
            {
                *(*o_buf) = '\0';
                return 0;
            }
            *(*o_buf)++ = min(count, 255);
            *(*o_buf)++ = code;
            *o_buflen -= 2;
            if (count <= 255)
            {
                return 1;
            }
            count -= 255;
        }
    }
    else
    {
        printf("invalid cigar format %d\n", format);
        return 1;
    }
}
int computeEditDistanceWithCigar(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, int useM,
    enum CigarFormat format)
{
    assert(k < MAX_K);

    int L[MAX_K+1][2 * MAX_K + 1];
    int i,j;
    for ( i = 0; i < MAX_K+1; i++)
    {
        for ( j = 0; j < 2*MAX_K+1; j++)
        {
            L[i][j] = -2;
        }
    }
    char A[MAX_K+1][2 * MAX_K + 1];

    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];

    const char* p = pattern;
    const char* t = text;
    if (NULL == text) return -1;            // This happens when we're trying to read past the end of the genome.

    int end = min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend)
    {
        uint64_t x = *((uint64_t*) p) ^ *((uint64_t*) t);
        if (x)
        {
            unsigned int zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end)
    {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
        if (useM)
        {
            if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M', format))
            {
                return -2;
            }
        }
        else
        {
            if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=', format))
            {
                return -2;
            }
            if (patternLen > end)
            {
                // Also need to write a bunch of X's past the end of the text
                if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X', format))
                {
                    return -2;
                }
            }
        }
        return 0;
    }
    int e = 0,d =0;
    for (e = 1; e <= k; e++)
    {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        for (d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d))
        {
            int best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best)
            {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best)
            {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d) + best;
            if (*p == *t)
            {
                int end = min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (1)
                {
                    uint64_t x = *((uint64_t*) p) ^ *((uint64_t*) t);
                    if (x)
                    {
                        unsigned int zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend)
                    {
                        best = end;
                        break;
                    }
                    t += 8;
                }
            }

            L[e][MAX_K+d] = best;

            if (best == patternLen)
            {
                // We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
                // trace back through the dynamic programming array to build up the CIGAR string.

                int straightMismatches = 0;
                for (i = 0; i < end; i++)
                {
                    if (pattern[i] != text[i])
                    {
                        straightMismatches++;
                    }
                }
                straightMismatches += patternLen - end;
                if (straightMismatches == e)
                {
                    // We can match with no indels; let's do that
                    if (useM)
                    {
                        //
                        // No inserts or deletes, and with useM equal and SNP look the same, so just
                        // emit a simple string.
                        //
                        if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M', format))
                        {
                            return -2;
                        }
                    }
                    else
                    {
                        int streakStart = 0;
                        int matching = (pattern[0] == text[0]);
                        for (i = 0; i < end; i++)
                        {
                            int newMatching = (pattern[i] == text[i]);
                            if (newMatching != matching)
                            {
                                if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'), format))
                                {
                                    return -2;
                                }
                                matching = newMatching;
                                streakStart = i;
                            }
                        }

                        // Write the last '=' or 'X' streak
                        if (patternLen > streakStart)
                        {
                            if (!matching)
                            {
                                // Write out X's all the way to patternLen
                                if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X', format))
                                {
                                    return -2;
                                }
                            }
                            else
                            {
                                // Write out some ='s and then possibly X's if pattern is longer than text
                                if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=', format))
                                {
                                    return -2;
                                }
                                if (patternLen > end)
                                {
                                    if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X', format))
                                    {
                                        return -2;
                                    }
                                }
                            }
                        }
                    }
                    return e;
                }

#ifdef TRACE_LV
                // Dump the contents of the various arrays
                printf("Done with e=%d, d=%d\n", e, d);
                for (int ee = 0; ee <= e; ee++)
                {
                    for (int dd = -e; dd <= e; dd++)
                    {
                        if (dd >= -ee && dd <= ee)
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++)
                {
                    for (int dd = -e; dd <= e; dd++)
                    {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
                int curD = d;
                int curE = e;
                for (curE = e; curE >= 1; curE--)
                {
                    backtraceAction[curE] = A[curE][MAX_K+curD];
                    if (backtraceAction[curE] == 'I')
                    {
                        backtraceD[curE] = curD + 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                    }
                    else if (backtraceAction[curE] == 'D')
                    {
                        backtraceD[curE] = curD - 1;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                    }
                    else     // backtraceAction[curE] == 'X'
                    {
                        backtraceD[curE] = curD;
                        backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                    }
                    curD = backtraceD[curE];
#ifdef TRACE_LV
                    printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD],
                           backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
                }

                int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
                if (useM)
                {
                    accumulatedMs = L[0][MAX_K+0];
                }
                else
                {
                    // Write out ='s for the first patch of exact matches that brought us to L[0][0]
                    if (L[0][MAX_K+0] > 0)
                    {
                        if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=', format))
                        {
                            return -2;
                        }
                    }
                }

                curE = 1;
                while (curE <= e)
                {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action)
                    {
                        actionCount++;
                        curE++;
                    }
                    if (useM)
                    {
                        if (action == '=' || action == 'X')
                        {
                            accumulatedMs += actionCount;
                        }
                        else
                        {
                            if (accumulatedMs != 0)
                            {
                                if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M', format))
                                {
                                    return -2;
                                }
                                accumulatedMs = 0;
                            }
                            if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action, format))
                            {
                                return -2;
                            }
                        }
                    }
                    else
                    {
                        if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action, format))
                        {
                            return -2;
                        }
                    }
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0)
                    {
                        if (useM)
                        {
                            accumulatedMs += backtraceMatched[curE];
                        }
                        else
                        {
                            if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=', format))
                            {
                                return -2;
                            }
                        }
                    }
                    curE++;
                }
                if (useM && accumulatedMs != 0)
                {
                    //
                    // Write out the trailing Ms.
                    //
                    if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M', format))
                    {
                        return -2;
                    }
                }
                *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
    return -1;
}
int landau2site(struct m_opt *opt,
                int length_r,int length_t,uint64_t r_start,
                 char *cigarBuf,unsigned int *siteBuf,
                 int mode,int scoreMode)
                 //mode 0  前向，mode 1 后向
{
    int start_cigar = 0;
    int t = 0,r = 0;
    int k = 0;
    int score = 0;

    r = r_start;
    k = 0;

    if(mode == 0)//前向
    {
        t = 0;
        while(cigarBuf[k]=='D')
        {
            start_cigar++;
            k++;
            r++;
        }
        if (scoreMode&&(start_cigar!=0)) score -= opt->miss;

        while(t<length_t)
        {
            if(cigarBuf[k]=='=')
            {
                score += opt->match;
                siteBuf[t]=r;
                t++;
                r++;
            }
            else if(cigarBuf[k]=='D')
            {
                if((cigarBuf[k-1]!='D')&&(k>0))
                    score -= opt->gap;
                r++;
            }
            else if(cigarBuf[k]=='I')
            {
                if((cigarBuf[k-1]!='I')&&(k>0))
                    score -= opt->miss;
                siteBuf[t]=r-1;
                t++;
            }
            else if(cigarBuf[k]=='X')
            {
                score -= opt->miss;
                siteBuf[t]=r;
                t++;
                r++;
            }
            k++;
        }
    }
    else if(mode == 1)
    {
        t = length_t-1;
        while(cigarBuf[k]=='D')
        {
            start_cigar++;
            k++;
            if(r!=0)r--;
        }
        if (scoreMode&&(start_cigar!=0)) score -= opt->miss;

        while(t>=0)
        {
            if(cigarBuf[k]=='=')
            {
                score += opt->match;
                siteBuf[t]=r;
                if(r!=0)r--;
                t--;
            }
            else if(cigarBuf[k]=='D')
            {
                if((cigarBuf[k+1]!='D')&&(k>0))
                    score -=opt->gap;
                if(r!=0)r--;
            }
            else if(cigarBuf[k]=='I')
            {
                if((cigarBuf[k+1]!='I')&&(k>0))
                    score -=opt->miss;
                siteBuf[t]=r;
                t--;
            }
            else if(cigarBuf[k]=='X')
            {
                score -= opt->miss;
                siteBuf[t]=r;
                if(r!=0)r--;
                t--;

            }
            k++;
        }
    }
    return score;
}

