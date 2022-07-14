import sys

seq_count = 0
base_name = '.'.join(sys.argv[1].split('.')[:-1])

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            if seq_count != 0:
                open_file.close()
            seq_count += 1
            open_file = open(base_name+'_'+str(seq_count)+'.fna', 'w')
            open_file.write(line)
        else:
            open_file.write(line)
    open_file.close()
