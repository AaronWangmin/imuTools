# from imuFactory import ImuFactory

class AllanVariance:
    def __init__(self,dataList,timeIntervel=1) -> None:
        self.dataList = dataList
        self.timeIntervel = timeIntervel      
        self.rms = 0

        self.caculateAllanVariance()
    
    def caculateAllanVariance(self):
        blockList = self.getSimpleBlock(self.dataList,self.timeIntervel)
        blockAverageList = self.getBlockAverageList(blockList)
        blockAverageDiffList = self.getDiffList(blockAverageList)        
        self.rms = self.getRms(blockAverageDiffList)
      
    # step 1: Blcok
    def getSimpleBlock(self,dataList,timeIntervel):
        blockList = []

        i = 0
        while i < len(dataList):       
            dataBlock = []
            for j in range(timeIntervel):                
                index = i + j 
                if(index < len(dataList)):
                    dataBlock.append(dataList[index])    

            i = i + timeIntervel
            blockList.append(dataBlock)

        return blockList         
    
    # step 2: average
    def getBlockAverageList(self,blockList):
        blockAverageList = []        
        for i in blockList:
            average = self.getDataListAverage(i)
            blockAverageList.append(average)  

        return blockAverageList     
    
    def getDataListAverage(self,dataList):
        sum = 0.0
        for i in dataList:
            sum += float(i)

        return sum/len(dataList)

    # step 3: difference
    def getDiffList(self,blockAverageList):
        blockAverageDiffList = []       
        i = 0
        while i < len(blockAverageList):
            if (i + 1) < len(blockAverageList):
                diff = blockAverageList[i] - blockAverageList[i+1]
                blockAverageDiffList.append(diff)            
            i += 1        
        return blockAverageDiffList

    # step 4: rms
    def getRms(self,blockAverageDiffList):
        average = self.getDataListAverage(blockAverageDiffList)
        
        averageDiffPow2List = []
        for i in blockAverageDiffList:
            averageDiffPow2List.append(pow(i - average,2))

        rms = pow(self.getDataListAverage(averageDiffPow2List),0.5)
        return rms


#--------------------------------test-------------     
# fd = ImuFactory()
# fd.getDataFromFile("E:\\project_201903\\ecarx-tools-202106\\imuProcess\\20210804150431_WANGMIN_SHANGHAI_AFA1119(1).ASC")

# xAcce = AllanVariance(fd.XAccelList,10000)
# print(xAcce.rms)
# print("hello")




                
        


        


