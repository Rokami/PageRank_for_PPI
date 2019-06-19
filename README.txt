1.主要程序为pagerank.py
2.各文件生成步骤在step.sh
3.流程和一些数据信息也在step.sh里
4.将所有文件放入一个文件夹才能运行成功
5.output.txt为epsilon=1e-7，alpha=0.2的beast cancer(OMIM ID=114480)结果
6.直接运行命令：python pagerank.py > <自己命名的文件名> 或者：python pagerank.py <接OMIM_ID> > <自己命名的文件名>.\
 如 python pagerank.py 211980 > LUNG_CANCER.txt