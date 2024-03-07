import os

vec=["cn0","cn1","cn2","cn3","cn4","cn5","cn6","cn7","cn8"]
for cn in vec:
    os.system("scp /root/pjt/GLEX_Coll_lib/build/test/* root@"+cn+":/root/")