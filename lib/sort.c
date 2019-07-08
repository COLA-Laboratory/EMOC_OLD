#include "../headers/global.h"
#include "../headers/sort.h"


extern void bublesort_weight(Weight_distance_info_t* distanceInfo, int size)
{
    int i = 0, j = 0;
    int temp_index = 0;
    double temp_distance;

    for(i=0;i<size;i++) //进行10次循环
    {
        for (j = i + 1; j < size; j++) //循环比较剩余的变量
        {
            if (distanceInfo[i].E_distance > distanceInfo[j].E_distance) //如果前面一个数比后面数大，交换两个数的值
            {
                temp_distance = distanceInfo[i].E_distance;
                temp_index = distanceInfo[i].idx;
                distanceInfo[i].E_distance = distanceInfo[j].E_distance;
                distanceInfo[i].idx = distanceInfo[j].idx;
                distanceInfo[j].idx = temp_index;
                distanceInfo[j].E_distance = temp_distance;
            }
        }
    }

}