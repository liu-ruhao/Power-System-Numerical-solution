readme
1，本代码目的在于 基于电力系统暂态稳定性分析的数值解法之分步计算法，
	计算极限切除角下的估计极限切除时间
2，给出必要参数后可以根据此代码求出估计的极限切除角，估计值与极限值的误差一般在3%以内
3，在cofig中输入数据，通过修改合理的 步长h 和计算时间 tm可以调整估计极限切除时间的精度
4，本代码并没有体现故障切除后的攻角与时间规律，仅考虑故障切除前的规律