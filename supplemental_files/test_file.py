with open('subject1.fasta', 'r') as file:
    for line in file:
        print(line.strip())  # 去除换行符并打印每一行内容
