class Triplets:
    def __init__(self):
        self.data = ([], ([], []))
        self.size_data = 0

    def __str__(self):
    	#revoir affichage
    	chaine="\tsize: "+str(self.size_data)+"\n"
    	for i in range (0,self.size_data):
    		chaine=chaine+"val:"+str(self.data[0][i])
    		chaine=chaine+" (I:"+str(self.data[1][0][i])+" J:"+str(self.data[1][1][i])+")\n"
    	return chaine


    def append(self, I, J, val):
    	self.size_data = self.size_data+1
    	self.data[0].append(val)
    	self.data[1][0].append(I)
    	self.data[1][1].append(J)
    	
