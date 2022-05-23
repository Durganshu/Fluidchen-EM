#include <iostream>
#include <string>

#include "Case.hpp"

void printIntro();

int main(int argn, char **args) {

    if (argn > 1) {
        printIntro();
        std::string file_name{args[1]};
        Case problem(file_name, argn, args);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}

void printIntro() {
    std::cout << "\n\n\n\nWelcome to FluidChen Flow Solver!\nThis project wqas developed as a part of Computational "
                 "Fluid Dynamics Lab course.\n\n\n\n";
    std::cout << "                                                     ..::::::..\n"
              << "                                                 .^!7?7^:::.:^^:\n"
              << "                                              .~?55YY!..:~:   .\n"
              << "                                            :7Y5GGP5YJ?77:...:^^^.\n"
              << "                                          :?Y5YG#BGG5J7~^::::^~~~:\n"
              << "                                        .7Y5YYYB&#P7:.\n"
              << "                                       :J55YYYYB#!\n"
              << "                                      ^Y5YYYYYYP~\n"
              << "                                     :Y5YYYYYYYJ.\n"
              << "                                     J5YYYYYYYJ?!..::::.\n"
              << "                                    ~YYYYYYYYYJ?7~^:::::\n"
              << "                                    7YYYYYYYYJJ?77^\n"
              << "                                    ?JJJYYYYJJJ???!\n"
              << "                                    ?JJJJJJJJJJJJJJ~\n"
              << "                                    !???JJJJJJYYYYYY^ \n"
              << "                                    ^????JJJYYYY555P5^\n"
              << "                                    .????JJJYY55PPPPGJ\n"
              << "                                     ^???JJY55PPPGGGBY\n"
              << "                                      !JJJYY5PPPGGGBP:\n"
              << "                                       ^?Y55PPPGGBGJ:\n"
              << "                                         :~7JJYJ?!:\n\n\n"
              <<

        "               ?G???7. 5J     JY    J5  YY .P5?Y57   :!~^^^~^  ?.   !~  ?!:^^:  Y7.   J\n"
              << "               JG      PY     Y5    YP  55 .B7  :GY .J.        J~:::?~  ?~.::   J:7~  J.\n"
              << "               JB777:  PY     YG.   P5  55 .B7  .PY .J.        J~:^:?~  ?!:::   J  ~7:J.\n"
              << "               JP      5P???! :Y5??Y5:  5Y .G5?J5?.  :7~^:^~^  J.   !~  ?~::::  J   :7Y\n\n\n\n";
}