# Bacterial-Gene-Clusters-Analysis-Tool

The 'ClusterFinding.py' file is the runnable file.

Instruction:
1. Open the termial, and set to the working directory, type in: python ClusterFinding.py
2. Then, the the program will require you to input name of the genome file you want to work on, the input should be a fatsta file of a bacterial's genome, the example input is the file: 'sequence.txt'
3. Then, input the protein sequence of the core gene you want to search for, the example protein sequence is: 
MAAARERAAGHRRRFDLRFEAYFDALPERLDTPALSRFTPRCLELLRDLSLRGGKRLRVALLYEAARLVTTGPVPGLAEAALSIELLQTHGLVHDDIIDDAPLRRGGPSTYYAYRQEFPAADRTALGLAVLAGDLAAFLSMRVLLEADVPAELRQAMLEVQLDAGAETVAGQIVDLERDLRRIPDEEFLHTVCEYKSTRYSVLAPLRLGLLAAGERLGGEEDARLRRYAIPVGIGGQLRDDYLDLFGDEDSTGKSTGADLRAGRRSYAVSALLAAADEEQRAVVESALGDPHCPAATVERVRELARRLGVDRKLRADMRRCAEAAAAEAGTWRPHWREEAVAFFAGLPLANVVGS
4. Then, input the percentage of the similarity of protein sequence comparison you want to set. The suggested similarity percentage is 30
5. Then, input the cluster length you want to set. The suggested length is 15000
6. After the running of the program, the the program will create a .txt file containing all the annotated genes in the found clusters names 'gene_cluster.txt'

