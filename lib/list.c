#include "../headers/list.h"


/*list*/
extern void list_del(struct list_head *prev, struct list_head *next)
{
    prev->next = next;
    next->prev = prev;
}

extern void list_del_init(struct list_head *head)
{
    list_del(head->prev, head->next);
    LIST_INIT_HEAD(head);
}


extern void list_add(struct list_head *new,
                     struct list_head *prev,
                     struct list_head *head)
{
    prev->next = new;
    head->prev = new;
    new->next = head;
    new->prev = prev;
}

extern void list_add_head(struct list_head *new, struct list_head *head)
{
    list_add(new, head, head->next);
}
extern void list_add_tail(struct list_head *new, struct list_head *head)
{
    list_add(new, head->prev, head);
}

extern void list_merge(struct list_head *head1, struct list_head *head2)
{
    head1->prev->next = head2->next;
    head2->next->prev = head1->prev;
    head1->prev = head2->prev;
    head2->prev->next = head1;

    LIST_INIT_HEAD(head2);
    return;
}


/*list*/