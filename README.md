# Joint Load and Capacity Scheduling for Flexible Radio Resource Management of High-Throughput Satellites

This work referred to	T. Ramirez, C. Mosquera and N. Alagha, Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads, IEEE Transactions on Broadcasting vol. 68, no. 3, pp. 723-739, 2022.
and their code: https://github.com/tomramzp/Resource.

First, familiarizing yourself with the aforementioned literature will facilitate your comprehension of this code.  The software requirements and main scripts are identical.

The following code allows replicating the results from the publication:
https://www.sciopen.com/article/10.26599/TST.2024.9010161

首先非常感谢开源代码的作者们，没有他们也没有今天这篇文章。

运行本文的程序需要安装软件：
“# SOFTWARE REQUIREMENTS

- MATLAB 2019b or newer version is required. The following MATLAB toolboxes are required:
	- Global Optimization toolbox
	- Parallel Computing Toolbox
	- Optimization Toolbox
	- Statistics and Machine Learning Toolbox

- Additionally, two external packages are also required, MOSEK and CVX: 

	- MOSEK ApS is a software for nonlinear convex optimization. The software requires a user license which can be obtained free of charge for researchers, students or other academic members. More information about the license can be found at https://www.mosek.com/products/academic-licenses/.  MOSEK is employed through CVX ( that is detailed below). It should be reminded to follow the steps in order as described in the guide. In particular,  Step 4  is required to detect the MOSEK licence with CVX, after placing the MOSEK license in the correct folder. 

	- CVX is a powerful and efficient solver for convex optimization problems. A user guide and installation steps can be found in http://web.cvxr.com/cvx/doc/mosek.html
       	  to install both MOSEK and CVX in MATLAB.
”
引用于https://github.com/tomramzp/Resource

各程序介绍：

System_Simulations_TwoDimesion_r7_paper1 主程序1
ResourceAssignmentTwoDimesion_r7_paper1 类与函数1
System_Simulations_TwoDimesion_r7_China 主程序2
ResourceAssignmentTwoDimesion_r7_China 类与函数2

draw_circle, get_circle_center,  get_distance_square, get_sign, min_cover_circle 画出一个覆盖给定点集中所有点的最小圆。
程序来自：https://blog.csdn.net/maple_2014/article/details/107973665

PlotResults_CR2 运行完主程序后，打开System simulations mm_dd_yyyy_mm_ss文件夹，可以看到数据已保存在场景文件夹中。运行PlotResults_CR2, 选择需要画图的场景文件夹即可。
当需要绘制 System_Simulations_TwoDimesion_r7_China 的结果时，将%% just for 64 beams部分的注释去除，%% for 6 beams部分加上注释即可。主程序中的绘图同理，如果需要画图，需要清除和绘图相关的所有注释。

