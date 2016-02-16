#include <stdlib.h>
#include <R.h>

int uniform_distribution(int rangeLow, int rangeHigh)
{
     double myRand = rand()/(1.0 + RAND_MAX);
     int range = rangeHigh - rangeLow + 1;
     int myRand_scaled = (myRand * range) + rangeLow;
     return myRand_scaled;
}

void draw_card(int *xc, int *x_ace)
{
     int c = uniform_distribution(1, 13);
     if (c > 10) c = 10;

     (*xc) += c;

     if (!(*x_ace) && c == 1)
     {
          (*xc) += 10;
          (*x_ace) = 1;
     }

     if ((*x_ace) && (*xc) > 21)
     {
          (*xc) -= 10;
          (*x_ace) = 0;
     }
}

void bj_transition(int *pc, int *p_ace, int *dc, int *d_ace, int *r, int a, int *sf)
{
     int stick = 1;
     int hit = 2;

     if (a == hit)
     {
          draw_card(pc, p_ace);

          if (*pc > 21)
          {
               *r = -1;
               *sf = 200; // state 200: agent loses
          }
          else *r = 0;
     }
     else if (a == stick)
     {
          // Dealer's fixed strategy: hit until 17+
          while (*dc < 17) draw_card(dc, d_ace);

          if (*dc > 21)
          {
               *r = 1;
               *sf = 202; // state 202: agent wins
          }
          else
          {
               int p_diff = 21 - *pc;
               int d_diff = 21 - *dc;

               if (p_diff < d_diff)
               {
                    *r = 1;
                    *sf = 202; // state 202: agent wins
               }
               else if (p_diff > d_diff)
               {
                    *r = -1;
                    *sf = 200; // state 200: agent loses
               }
               else
               {
                    *r = 0;
                    *sf = 201; // state 201: draw
               }
          }
     }

}

unsigned int state_index(int pc, int p_ace, int dc)
{
     return (p_ace * 100) + ((pc - 12) * 10) + (dc - 2);
}

void bj_evaluate_policy(int *pi, int *num_episodes, double *ret)
{

     *ret = 0.0;

     int i = 0;
     for (i = 0; i < *num_episodes; ++i)
     {
          int pc = 0;
          int p_ace = 0;
          int dc = 0;
          int d_ace = 0;
          int sf = 0;

          while (pc < 12) draw_card(&pc, &p_ace);

          draw_card(&dc, &d_ace);

          int r = 0;
          while (sf == 0)
          {
               int a = pi[state_index(pc, p_ace, dc)];
               bj_transition(&pc, &p_ace, &dc, &d_ace, &r, a, &sf);
          }

          (*ret) += r;
     }

     (*ret) /= (float)(*num_episodes);
}
