#include "util.h"
#include <ctime>

/*
* Generate random number that represents a card (1-13)
*/
int generate_card()
{
  return rand() % 13 + 1;
}

/*
* The AI will always have a card between 15-21. This will generate their card.
*/
int generate_ai_card()
{
  return rand() % 7 + 15;
}

int main()
{
  return 0;
}
