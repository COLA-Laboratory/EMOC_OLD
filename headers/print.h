# ifndef _PRINT_H_
# define _PRINT_H_

# include "../headers/global.h"
extern void print_variable (char *file_name, SMRT_individual *pop_table);
extern void print_objective (char *file_name, SMRT_individual * pop_table);
extern void print_error (int condition, int n, ...);
extern void print_progress ();
# endif