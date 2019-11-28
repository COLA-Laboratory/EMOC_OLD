#ifndef SAMARITAN_VECTOR_H
#define SAMARITAN_VECTOR_H


typedef struct int_vector
{
    int value;
    struct int_vector* next;
}i_vector;


typedef struct double_vector
{
    double value;
    struct double_vector* next;
}d_vector;

int int_vector_size (struct int_vector* head);
void int_vector_pushback (struct int_vector* head, int value);
int int_vector_pop (struct int_vector* head);
int int_vector_get (struct int_vector* head, int index);
void int_vector_set (struct int_vector* head, int index, int value);
void int_vector_free (struct int_vector* head);
void int_vector_print (struct int_vector* head);
void int_vector_remove(struct int_vector* head, int index);

int double_vector_size (struct double_vector* head);
void double_vector_pushback (struct double_vector* head, double value);
double double_vector_pop (struct double_vector* head);
double double_vector_get (struct double_vector* head, int index);
void double_vector_free (struct double_vector* head);
void double_vector_print (struct double_vector* head);

#endif //SAMARITAN_VECTOR_H