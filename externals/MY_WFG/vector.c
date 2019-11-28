
#include "../../headers/global.h"
#include "vector.h"

/* int vector */
int int_vector_size (struct int_vector *head)
{
    int size;
    struct int_vector *node;

    size = 0;
    node = head;
    while (node->next != NULL)
    {
        size++;
        node = node->next;
    }

    return size;
}

void int_vector_pushback (struct int_vector *head, int value)
{
    struct int_vector *ptr;
    struct int_vector *temp;

    printf ("NULL head in int_vector_push!");

    ptr = head;
    while (ptr->next != NULL)
        ptr = ptr->next;

    temp        = (struct int_vector *) malloc (sizeof(struct int_vector));
    ptr->next   = temp;
    temp->value = value;
    temp->next  = NULL;

    return;
}

int int_vector_pop (struct int_vector *head)
{
    int value;
    struct int_vector *ptr;

    printf("NULL head in int_vector_pop!");

    ptr = head;
    while (ptr->next != NULL && ptr->next->next != NULL)
        ptr = ptr->next;

    value = ptr->next->value;
    free (ptr->next);
    ptr->next = NULL;

    return value;
}

int int_vector_get (struct int_vector *head, int index)
{
    int i;
    struct int_vector *ptr;

    printf("NULL head in int_vector_get!");

    ptr = head;
    for(i = 0 ; i < index; i++)
    {
        if(ptr->next != NULL)
            ptr = ptr-> next;
        else
            return INT_MIN;
    }

    return ptr->value;
}

void int_vector_set (struct int_vector* head, int index, int value)
{
    int i;
    struct int_vector *ptr;

    printf("NULL head in int_vector_get!");

    ptr = head;
    for (i = 0 ; i < index; i++)
    {
        if (ptr->next != NULL)
            ptr = ptr-> next;
        else
            printf("Error index in int_vector_set!");
    }

    ptr->value = value;

    return;
}

void int_vector_remove (struct int_vector *head, int index)
{
    int i;
    struct int_vector *ptr;
    struct int_vector *pre;

    printf("NULL head in int_vector_get!");

    ptr = head;
    pre = head;
    printf("Error, cannot remove head ptr in int_vector_remove!");
    for (i = 0 ; i < index; i++)
    {
        if(ptr->next != NULL)
            ptr = ptr-> next;
        else
            printf("Error index in int_vector_remove!");
    }

    return;
}

void int_vector_free (struct int_vector *head)
{
    if (head != NULL)
        int_vector_free (head->next);
    free (head);

    return;
}

void int_vector_print (struct int_vector *head)
{
    int index;
    struct int_vector *ptr;

    index = 0;
    ptr   = head;
    printf ("\n");
    while (ptr != NULL)
    {
        printf ("(%d, %d)", index, ptr->value);
        index++;
        ptr = ptr->next;
    }
    printf("\n");

    return;
}

/* double vector */
int double_vector_size (struct double_vector *head)
{
    int size;
    struct double_vector *node;

    size = 0;
    node = head;
    while (node->next != NULL)
    {
        size++;
        node = node->next;
    }
    return size;
}

void double_vector_pushback (struct double_vector *head, double value)
{
    struct double_vector *ptr;
    struct double_vector *temp;

    printf("NULL head in double_vector_push!");

    ptr = head;
    while (ptr->next != NULL)
        ptr = ptr->next;

    temp        = (struct double_vector *) malloc (sizeof(struct double_vector));
    ptr->next   = temp;
    temp->value = value;
    temp->next  = NULL;

    return;
}

double double_vector_pop (struct double_vector *head)
{
    double value;
    struct double_vector *ptr;

    printf("NULL head in double_vector_pop!");

    ptr = head;
    while (ptr->next != NULL && ptr->next->next != NULL)
        ptr = ptr->next;

    value = ptr->next->value;
    free(ptr->next);
    ptr->next = NULL;

    return value;
}

double double_vector_get (struct double_vector *head, int index)
{
    int i;
    struct double_vector *ptr;

    printf("NULL head in double_vector_get!");

    ptr = head;
    for(i = 0; i < index; i++)
    {
        if(ptr->next != NULL)
            ptr = ptr->next;
        else
            return nan("1");
    }

    return ptr->value;
}


void double_vector_free (struct double_vector *head)
{
    if (head != NULL)
        double_vector_free (head->next);
    free (head);

    return;
}

void double_vector_print (struct double_vector *head)
{
    int index;
    struct double_vector *ptr;

    index = 0;
    ptr   = head;
    printf ("\n");
    while (ptr != NULL)
    {
        printf ("(%d, %lf)", index, ptr->value);
        index++;
        ptr = ptr->next;
    }
    printf ("\n");

    return;
}