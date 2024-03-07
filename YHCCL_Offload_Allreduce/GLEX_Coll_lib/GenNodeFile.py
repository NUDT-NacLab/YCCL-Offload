import commands 
import os
import sys

def fetch_idle_node_list(PartitionN,Rack_ID):
    env_dist = os.environ
    # status,output = commands.getstatusoutput("yhinfo -p "+PartitionN+" --state=idle")
    status,output = commands.getstatusoutput("yhinfo -M " + Rack_ID+ " -p " + PartitionN + " --state=idle" )
    start=output.find('[')
    end=output.find(']') 
    tmp = (output[start+1:end]).split(',')
    #print(tmp)
    IdleNodeVec=[]
    for nodestr in tmp:
        t = nodestr.find('-')
        if t==-1:
            IdleNodeVec.append(int(nodestr))
        else:
            start=int(nodestr[:t])
            end  =int(nodestr[t+1:])
            for nodeid in range(start,end+1):
                IdleNodeVec.append(nodeid)
    return IdleNodeVec
def main(argv):
    N = int(argv[1])
    if(len(argv) >= 4):
        Rack_ID=argv[3]
        print(Rack_ID)
    PartitionN=argv[2]
    print(PartitionN)
    Contiguous = 1
    IdleLists = fetch_idle_node_list(PartitionN,Rack_ID)
    step = len(IdleLists)/N
    fileN = "nodefile.txt"
    f = open(fileN,'w')
    restr = ""
    i = 0
    while i < N:
        step = min(Contiguous,N-i)
        for j in range(0,step):
            if i+j != 0:
                restr += (",cn"+str(IdleLists[(i*step+j)]))
            else:
                restr += ("cn"+str(IdleLists[(i*step+j)]))
        i+=step
    f.write(restr)
    f.close()
 
if __name__ == '__main__':
    main(sys.argv)