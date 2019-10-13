#ifndef _EA_LIB_H_
#define _EA_LIB_H_


#define LIST_INIT_HEAD(ptr)  \
    (ptr)->next = (ptr); (ptr)->prev = (ptr);


#define list_for_each_safe(pose, safe, head) \
    for((pose) = (head)->next, (safe) = (pose)->next; (pose) != (head); (pose) = (safe), (safe) = (pose)->next)


#define list_for_each(pose, head) \
    for ((pose) = (head)->next; (pose) != (head); (pose) = (pose)->next)

#define list_empty(ptr) \
    (ptr)->next == (ptr)

struct list_head{
    struct list_head *prev, *next;
};
extern void list_del(struct list_head *prev, struct list_head *next);
extern void list_del_init(struct list_head *head);
extern void list_add(struct list_head *new, struct list_head *prev, struct list_head *head);
extern void list_add_head(struct list_head *new, struct list_head *head);
extern void list_add_tail(struct list_head *new, struct list_head *head);
extern void list_merge(struct list_head *head1, struct list_head *head2);



#endif