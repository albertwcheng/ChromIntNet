#!/bin/env python

from sys import *
from operator import attrgetter
import os.path
from getopt import getopt


class SuperCluster:
    def __init__(self,_chrom,_start0,_end1):
        self.chrom=_chrom
        self.start0=_start0
        self.end1=_end1
        self.memberClusters=dict() #[cell]->[nodeClusters]
        self.name=""

    def addCluster(self,cluster):
        cellName=cluster.getCellName()
        try:
            memberClusterList=self.memberClusters[cellName]
        except:
            memberClusterList=[]
            self.memberClusters[cellName]=memberClusterList

        memberClusterList.append(cluster)
        cluster.superClusterParent=self

class MultiChromNet:
    def __init__(self):
        self.superClusters=dict() #[chrom]->[SuperCluster]
        self.chromIntNets=dict() #[cell]->ChromIntNet
        self.superClusterPrefix="SuperCluster_"
        self.superClusterIndex=0
        self.allClusters=dict() #[chrom]->[cluster]  #need to be sorted after loading every file
        self.tempClusterList=[]
        #self.interactionType="ChIAPETLoop"

    def loadChIAPETFile(self,cellName,filename):
        cinet=ChromIntNet()
        cinet.cellName=cellName
        cinet.clusterPrefix=cellName+"_"
        cinet.readChIAPETFile(filename)
        cinet.finalizeNodeAndNodeClusters()
        self.chromIntNets[cellName]=cinet
        
        for chrom,clusterList in cinet.nodeClusters.items():
            try:
                allClusterList=self.allClusters[chrom]
            except:
                allClusterList=[]
                self.allClusters[chrom]=allClusterList

            allClusterList.extend(clusterList)

        #now form supercluster for individual clusters

    def formSuperCluster(self,superClusterList,chrom,clusterBoundStart0,clusterBoundEnd1):
        if len(self.tempClusterList)==0:
            return

        self.superClusterIndex+=1
        
        supercluster=SuperCluster(chrom,clusterBoundStart0,clusterBoundEnd1)
        supercluster.name=self.superClusterPrefix+str(self.superClusterIndex)
    
        superClusterList.append(supercluster)
        
        for cluster in self.tempClusterList:
            supercluster.addCluster(cluster)
    
        self.tempClusterList=[]
        

    def finalizeSuperClusters(self):
        self.tempClusterList=[]
        
        #first sort cluster
        #form clusters noting the cellName
        for chrom,allClusterList in self.allClusters.items():
            allClusterList.sort(key=attrgetter('start0'))
        
            superClusterList=[]
            self.superClusters[chrom]=superClusterList
        
            clusterBoundStart0=allClusterList[0].start0
            clusterBoundEnd1=allClusterList[1].end1
            for i in range(0,len(allClusterList)):
                thisCluster=allClusterList[i]
                if thisCluster.start0>=clusterBoundEnd1:
                    self.formSuperCluster(superClusterList,chrom,clusterBoundStart0,clusterBoundEnd1)
                    clusterBoundStart0=thisCluster.start0
                    clusterBoundEnd1=thisCluster.end1
    
                clusterBoundEnd1=max(clusterBoundEnd1,thisCluster.end1)
                self.tempClusterList.append(thisCluster)
            
            self.formSuperCluster(superClusterList,chrom,clusterBoundStart0,clusterBoundEnd1)
    

    def printSummary(self,ostream):
        print >> ostream,"ChromIntNets:"

        for cellName,cinet in self.chromIntNets.items():
            print >> ostream,"ChromIntNet "+cellName
            cinet.printSummary(ostream)
            print >> ostream,"---------"
        
        print >> ostream,"AllClusters:"

        for chrom,allClusterList in self.allClusters.items():
            for cluster in allClusterList:
                print >> ostream,cluster.getCellName()+"@"+cluster.name+" "+cluster.chrom+":"+str(cluster.start0+1)+"-"+str(cluster.end1)



        print >> ostream,"SuperClusters:"
        for chrom,superClusterList in self.superClusters.items():
            for supercluster in superClusterList:
                print >> ostream, supercluster.name+" "+supercluster.chrom+":"+str(supercluster.start0+1)+"-"+str(supercluster.end1)
                for cellName,memberClusterList in supercluster.memberClusters.items():
                    print >> ostream,"\t"+cellName+" :",
                    for memberCluster in memberClusterList:
                        print >> ostream,memberCluster.name+" "+memberCluster.chrom+":"+str(memberCluster.start0+1)+"-"+str(memberCluster.end1),"[",
                        for node in memberCluster.nodes:
                            print >> ostream,str(node.start0+1)+"-"+str(node.end1),
                        print >> ostream,"]",
                    print >> ostream,""

    def writeClusterNet(self,folder):
        filClusterNetSif=open(folder+"clusterNet.sif","w")
        filClusterNetEdgeAttributeTable=open(folder+"clusterNet.edgeAttrs.txt","w")
        filClusterBed=open(folder+"superclusters.bed","w")
        filNodesBed=open(folder+"nodes.bed","w")
        
        print >> filClusterNetEdgeAttributeTable,"name\tInteractionStrength\tInteractionDegree"
        
        
        for chrom,superClusterList in self.superClusters.items():
            for supercluster in superClusterList:
                #output clusterbed
                print >> filClusterBed,supercluster.chrom+"\t"+str(supercluster.start0)+"\t"+str(supercluster.end1)+"\t"+supercluster.name   #+"\t"+str(len(cluster.nodes))
                
                strengthDict=dict() #[SuperClusterName]->[Cell]->Strength
                degreeDict=dict()  #[SuperClusterName]->[Cell]->Degree
                
                for cellName,clusterList in supercluster.memberClusters.items():
                    for cluster in clusterList:
                        for node in cluster.nodes:
                            
                            #output nodebed
                            print >> filNodesBed,node.chrom+"\t"+str(node.start0)+"\t"+str(node.end1)+"\t"+supercluster.name+"/"+cellName+"/"+cluster.name+"/"+node.name

                            if node.next:

                                connectingCluster=node.next.parent
                                connectingSuperCluster=node.next.parent.superClusterParent
                                connectingSuperClusterName=connectingSuperCluster.name
                                
                                try:
                                    strengthDictDict=strengthDict[connectingSuperClusterName]
                                    degreeDictDict=degreeDict[connectingSuperClusterName]
                                except:
                                    strengthDictDict=dict()
                                    degreeDictDict=dict()
                                    strengthDict[connectingSuperClusterName]=strengthDictDict
                                    degreeDict[connectingSuperClusterName]=degreeDictDict
                                
                                
                                if cellName not in strengthDictDict:
                                    strengthDictDict[cellName]=0.0
                                    degreeDictDict[cellName]=0
                                
                                
                                strengthDictDict[cellName]+=node.strength
                                degreeDictDict[cellName]+=1
        
                #print >> stderr,strengthDict
                
                #now output ClusterNet files
                for connectingSuperClusterName,strengthDictDict in strengthDict.items():
                    degreeDictDict=degreeDict[connectingSuperClusterName]
                    for cellName,strength in strengthDictDict.items():
                        degree=degreeDictDict[cellName]
                        print >> filClusterNetSif,supercluster.name+" "+cellName+" "+connectingSuperClusterName
                        print >> filClusterNetEdgeAttributeTable,supercluster.name+" ("+cellName+") "+connectingSuperClusterName+"\t"+str(strength)+"\t"+str(degree)

        filClusterNetSif.close()
        filClusterNetEdgeAttributeTable.close()
        filClusterBed.close()


class BoundNode:
    def __init__(self,_name,_chrom,_start0,_end1):
        self.name=_name
        self.next=None
        self.strength=0.0
        self.chrom=_chrom
        self.start0=_start0
        self.end1=_end1
        self.parent=None

    def getSuperClusterParent(self):
        if self.parent:
            return self.parent.superClusterParent

        return None

    def biggerThan(self,anotherNode):
        if self.chrom<anotherNode.chrom:
            return False
        if self.chrom>anotherNode.chrom:
            return True
        #now this.chrom==anotherNode.chrom
        if self.start0==anotherNode.start0:
            return self.end1>anotherNode.end1

        return self.start0>anotherNode.start0

class BoundNodeCluster:
    def __init__(self,_name,_chrom,_start0,_end1,_nodes,_cinetParent):
        self.name=_name
        self.nodes=_nodes
        self.chrom=_chrom
        self.start0=_start0
        self.end1=_end1
        self.superClusterParent=None
        self.parentNet=_cinetParent
        for node in self.nodes:
            node.parent=self

    def getCellName(self):
        return self.parentNet.cellName



class ChromIntNet:
    def __init__(self):
        self.nodes=dict() #chrom->list
        self.nodeClusters=dict() #chrom->list
        self.tempNodeList=[]
        self.clusterPrefix="cluster_"
        self.clusterIndex=0
        self.interactionType="ChIAPETLoop"
        self.multinetParent=None
        self.CellName=""
    

    
    def printSummary(self,ostream):
        for chrom,clusterList in self.nodeClusters.items():
            for cluster in clusterList:
                print >> ostream,cluster.name+" "+cluster.chrom+":"+str(cluster.start0+1)+"-"+str(cluster.end1)
                for node in cluster.nodes:
                    print >> ostream,node.name+" "+node.chrom+":"+str(node.start0+1)+"-"+str(node.end1),
                    if node.next:
                        print >> ostream,"("+str(node.strength)+") "+node.next.parent.name+"/"+node.next.name+" "+node.next.chrom+":"+str(node.next.start0+1)+"-"+str(node.next.end1)
                    else:
                        print >> ostream,""

    
    
    def formCluster(self,clusterList,chrom,clusterBoundStart0,clusterBoundEnd1):
        if len(self.tempNodeList)==0:
            return
        
        self.clusterIndex+=1
        
        cluster=BoundNodeCluster(self.clusterPrefix+str(self.clusterIndex),chrom,clusterBoundStart0,clusterBoundEnd1,self.tempNodeList,self)
        
        #try:
        #    clusterList=self.nodeClusters[chrom]
        #except:
        #    clusterList=[]
        #    self.nodeClusters[chrom]=clusterList
        
        clusterList.append(cluster)
        
        self.tempNodeList=[]
    
    

    def finalizeNodeAndNodeClusters(self):
        #sort the node
        
        self.tempNodeList=[]
        
        for chrom,nodeList in self.nodes.items():
            nodeList.sort(key=attrgetter('start0'))
            
            clusterList=[]
            self.nodeClusters[chrom]=clusterList
            
            #now form cluster for this chromosome
            clusterBoundStart0=nodeList[0].start0
            clusterBoundEnd1=nodeList[0].end1
            for i in range(0,len(nodeList)):
                thisNode=nodeList[i]
                if thisNode.start0>=clusterBoundEnd1:
                    #not continued
                    self.formCluster(clusterList,chrom,clusterBoundStart0,clusterBoundEnd1)
                    clusterBoundStart0=thisNode.start0
                    clusterBoundEnd1=thisNode.end1

                clusterBoundEnd1=max(clusterBoundEnd1,thisNode.end1)
                self.tempNodeList.append(thisNode)
            
            #end of the whole chromosome,form cluster with the last
            self.formCluster(clusterList,chrom,clusterBoundStart0,clusterBoundEnd1)


    def readChIAPETFile(self,filename):
        fil=open(filename)
        for lin in fil:
            lin=lin.rstrip()
            fields=lin.split("\t")
            name,chrom1,start1,end1,chrom2,start2,end2,strength=fields
            start1=int(start1)
            end1=int(end1)
            start2=int(start2)
            end2=int(end2)
            strength=float(strength)
            thisNode1=BoundNode(name+".1",chrom1,start1,end1)
            thisNode2=BoundNode(name+".2",chrom2,start2,end2)
            if thisNode1.biggerThan(thisNode2):
                thisNode2.next=thisNode1
                thisNode2.strength=strength
            else:
                thisNode1.next=thisNode2
                thisNode1.strength=strength

            try:
                nodeList=self.nodes[chrom1]
            except:
                nodeList=[]
                self.nodes[chrom1]=nodeList
            
            nodeList.append(thisNode1)
                
                
            try:
                nodeList=self.nodes[chrom2]
            except:
                nodeList=[]
                self.nodes[chrom2]=nodeList
            
            nodeList.append(thisNode2)
                
        
        fil.close()

    def writeClusterNet(self,folder):
        filClusterNetSif=open(folder+"clusterNet.sif","w")
        filClusterNetEdgeAttributeTable=open(folder+"clusterNet.edgeAttrs.txt","w")
        filClusterBed=open(folder+"clusters.bed","w")
        filNodesBed=open(folder+"nodes.bed","w")
        
        print >> filClusterNetEdgeAttributeTable,"name\tInteractionStrength\tInteractionDegree"
        

        for chrom,clusterList in self.nodeClusters.items():
            for cluster in clusterList:
                #output clusterbed
                print >> filClusterBed,cluster.chrom+"\t"+str(cluster.start0)+"\t"+str(cluster.end1)+"\t"+cluster.name+"\t"+str(len(cluster.nodes))
                strengthDict=dict()
                degreeDict=dict()
                for node in cluster.nodes:
                    #output nodebed
                    print >> filNodesBed,node.chrom+"\t"+str(node.start0)+"\t"+str(node.end1)+"\t"+cluster.name+"/"+node.name
                    if node.next:
                        connectingCluster=node.next.parent
                        try:
                            strengthDict[connectingCluster.name]+=node.strength
                            degreeDict[connectingCluster.name]+=1
                        except:
                            strengthDict[connectingCluster.name]=node.strength
                            degreeDict[connectingCluster.name]=1
                    
                #now output ClusterNet files
                for connectingClusterName,strength in strengthDict.items():
                    print >> filClusterNetSif,cluster.name+" "+self.interactionType+" "+connectingClusterName
                    print >> filClusterNetEdgeAttributeTable,cluster.name+" ("+self.interactionType+") "+connectingClusterName+"\t"+str(strength)+"\t"+str(degreeDict[connectingClusterName])

        filClusterNetSif.close()
        filClusterNetEdgeAttributeTable.close()
        filClusterBed.close()

def printUsageAndExit(programName):
    print >> stderr,programName,"subProgram [opts] args"
    print >> stderr,""
    
    print >> stderr,programName,"ChIAPET2Net [options] interactionFile outputPrefix"
    print >> stderr,"--cluster-prefix prefix. Set prefix to the cluster names [cluster_]"
    print >> stderr,""
    
    print >> stderr,programName,"Multinet --output-prefix outputPrefix interactionFile1 interactionFile2 ... interactionFileN"
    print >> stderr,"--output-prefix prefix. Prefix of output files [mandatory option]"
    
    exit(0)


def ChIAPET2Net(programName,subProgramName,opts,args):
    cinet=ChromIntNet()
    for o,v in opts:
        if o=='--cluster-prefix':
            cinet.clusterPrefix=v

    try:
        interactionFile,folder=args
    except:
        printUsageAndExit(programName)


    cinet.readChIAPETFile(interactionFile)
    cinet.finalizeNodeAndNodeClusters()
    #cinet.printSummary(stdout)
    cinet.writeClusterNet(folder)


def Multinet(programName,subProgramName,opts,args):
    mnet=MultiChromNet()
    outputPrefix=None
    
    for o,v in opts:
        if o=='--output-prefix':
            outputPrefix=v

    if not outputPrefix:
        print >> stderr,"no outputPrefix defined. Abort"
        printUsageAndExit(programName)

    for arg in args:
        cellName,filename=arg.split("=")
        mnet.loadChIAPETFile(cellName,filename)

    mnet.finalizeSuperClusters()

    #mnet.printSummary(stdout)

    mnet.writeClusterNet(outputPrefix)


if __name__=='__main__':
    programName=argv[0]
    subProgramList={"ChIAPET2Net":ChIAPET2Net,"Multinet":Multinet}
    
    try:
        subProgramName=argv[1]
        opts,args=getopt(argv[2:],'',['cluster-prefix=','output-prefix='])
   
   

        if subProgramName not in subProgramList.keys():
            print >> stderr,"subprogram %s not found. Abort" %(subProgramName)
            print >> stderr,"program list",",".join(subProgramList.keys())
            printUsageAndExit(programName)
    except:
        printUsageAndExit(programName)




    subProgramList[subProgramName](programName,subProgramName,opts,args)


