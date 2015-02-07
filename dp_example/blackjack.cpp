#include <iostream>
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

/*
 * Let's play the game!
 */
int main()
{
  // This creates a random seed for rand() based on the time.

  srand(time(NULL));

  bool play;
  int ai_card, my_card;
  char x;

  std::cout << "Welcome to a very simple BlackJack game to kill time.\n";

  play = true;

  // Get the AI's card.
  ai_card = generate_ai_card();

  // Get your first card.
  my_card = generate_card();

  std::cout << "You have " << my_card << " - (h)it or (s)tay?: ";

  // Ask for a hit or stay.
  do
    {
      std::cin >> x;

      // If they have hit...
      if ( x == 'h' )
        {
          // Generate a card, and add it on.
          my_card = my_card + generate_card();

          // Have they bust?
          if ( my_card > 21 )
            {
              std::cout << "You got bust (" << my_card << ").\n";
              play = false;
            }
          else
            {
              std::cout << "You now have " << my_card << " - (h)it or (s)tay?: ";
            }
        }
      // If they stayed...
      else if ( x == 's' )
        {

          // ... See who is bigger.
          if ( my_card > ai_card )
            {
              // You WIN!
              std::cout << "You won! (" << my_card << " vs " << ai_card << ")\n";
              play = false;
            }
          else
            {
              // You have FAILED.
              std::cout << "You lost! (" << my_card << " vs " << ai_card << ")\n";
              play = false;
            }
        }
      else
        {
          std::cout << "Please either type \"h\" to hit or \"s\" to stay: ";
        }
    }
  while ( play == true );
  // Some credits and stuff.
  std::cout << "Thank you for playing Dion Moult's C++ BlackJack game.\n";
  return 0;
}
