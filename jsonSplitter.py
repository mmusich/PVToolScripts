import json

shortIOVlist=["290543", "294034", "297281", "297598", "298653", "299443", "300155", "300389", "300538", "301046", "301417", "302131", "303719", "303790", "303998", "304661", "304911", "305040", "305204", "305898", "306029"]
longIOVlist =["290543", "294034", "296702", "296966", "297224", "297281", "297429", "297467", "297484", "297494", "297503", "297557", "297598", "297620", "297660", "297670", "298653", "298678", "298996", "299062", "299096", "299184", "299327", "299368", "299381", "299443", "299480", "299592", "299594", "299649", "300087", "300155", "300233", "300237", "300280", "300364", "300389", "300399", "300459", "300497", "300515", "300538", "300551", "300574", "300636", "300673", "300780", "300806", "300812", "301046", "301417", "302131", "302573", "302635", "303719", "303790", "303825", "303998", "304170", "304505", "304661", "304672", "304911", "305040", "305113", "305178", "305188", "305204", "305809", "305842", "305898", "305967", "306029", "306042", "306126", "306169", "306417", "306459", "306460", "306705", "306826"]

input_file = open ('json_DCSONLY.txt')
d = json.load(input_file)
    #print k,v
    
for i in range(0,len(shortIOVlist)):
    output={}
    if(i<len(shortIOVlist)-1):
        print i,shortIOVlist[i],shortIOVlist[i+1]
        for k,v in d.items():
            if((int(k)>=int(shortIOVlist[i])) and (int(k)<int(shortIOVlist[i+1]))):
                output[k]=v
    else:
        print i,shortIOVlist[i],999999
        for k,v in d.items():
            if((int(k)>=int(shortIOVlist[i])) and (int(k)<999999)):
                output[k]=v

    print output
    with open('data_IOV'+shortIOVlist[i]+'.json', 'w') as outfile:
        json.dump(output, outfile)

#print json_array
