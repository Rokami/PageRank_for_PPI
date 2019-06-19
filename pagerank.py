# -*- coding: utf-8 -*-

from pygraph import UndirectedGraph
import sys
class PageRank(object):
    """docstring for PageRank"""

    """define variables"""
    def __init__(self, graph=UndirectedGraph(),s=[]):
        self.graph = graph
        self.epsilon=1e-8#####epsilon:can change
        self.alpha=0.2########alpha:can change
        self.r=s
        self.query=[]
        self.record={i:0 for i in self.graph.get_all_node_ids()}
        self.p={i:0 for i in self.graph.get_all_node_ids()}
        for i in self.r:
            if i['data']>=self.epsilon*len(i['edges']):
                self.query.append(i)
        # print len(self.query)
    """Push method"""
    def Push(self,u={}):
        temp_p=self.p[u['id']]+self.alpha*u['data']
        temp_data=(1-self.alpha)*u['data']/2
        # print len(graph.neighbors(u['id']))
        for vertex in graph.neighbors(u['id']):
            temp_neighbor=self.graph.get_node(vertex)
            temp_neighbor['data']=temp_neighbor['data']+(1-self.alpha)*u['data']/(2*len(u['edges']))
            if temp_neighbor['data']>=self.epsilon*len(temp_neighbor['edges']) and self.record[temp_neighbor['id']]==0:
                self.query.append(temp_neighbor)
                self.record[temp_neighbor['id']]=1
            self.graph.get_node(vertex)['data']=temp_neighbor['data']
        self.p[u['id']]=temp_p
        u['data']=temp_data
        return u
    """ApproximatePR method"""
    def ApproximatePR(self):
        num=0
        while self.query:
            num+=1
            node=self.query.pop(0)
            self.record[node['id']]=0
            # print len(self.query)
            # print node['data']
            node=self.Push(node)
            # print node['data']
            # print "\n"
            if node['data']>=self.epsilon*len(node['edges']) and self.record[node['id']]==0:
                self.query.append(node)
                self.record[node['id']]=1
            # print "iteration: ",num,"query length: ",len(self.query)
        print "iteration: ",num,"query length: ",len(self.query)
        return self.p

if __name__ == '__main__':
    graph=UndirectedGraph()
    node={}
    node_map={}
    reverse_node_map={}
    edge=[]
    omim={}
    disease=[]
    posibility={}
    s=[]
    protein_name={}
    try:
        query_id=sys.argv[1]
    except Exception as e:
        query_id='114480'
    ############################################################
    ###build graph
    with open('protein_protein_interactions.txt','r') as f:
        list1=f.readlines()
        for line in list1:
            info=line.strip().split()
            node[info[0]]=''
            # print info[0]
            # print info[1]
            node[info[1]]=''
            edge.append([info[0],info[1]])
    for i in sorted(node.keys()):
        node_map[i]=graph.next_node_id
        reverse_node_map[graph.next_node_id]=i
        new=graph.new_node()
        graph.get_node(new)['data']=0
    # print range(graph.generate_node_id())
    ###add graph edge
    for x in edge:
        graph.new_edge(node_map[x[0]],node_map[x[1]])
    ############################################################
    ###read diease-gene file
    with open('omim_HPRD.txt','r') as f:
        list2=f.readlines()
        for line in list2:
            info=line.strip().split('\t')
            if(len(info)==1):
                omim[info[0]]=''
            else:
                gene=info[1].split(',')
                omim[info[0]]=gene
            disease.append(info[0])
    ###read score of disease similarity
    with open('MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.mat','r') as f:
        list3=f.readlines()
        for line in list3:
            # print line,
            info=line.strip().split('\t')
            if info[0]==query_id:
                posibility={disease[i]:float(info[i+1]) for i in range(len(disease))}
                break
    ###read HPRD_ID to protein name file and get dict variables
    with open('all_protein_HPRD_ID.txt','r')as f:
        list4=f.readlines()
        for line in list4:
            info=line.strip().split('\t')
            protein_name[info[1]]=info[0]
    ########################################################################
    #get inital vector s
    sumup=float(0)
    for x in posibility.keys():
        if posibility[x]>=0.6:
            if len(omim[x])==0:
                next
            for name in omim[x]:
                sumup=sumup+posibility[x]
                graph.get_node(node_map[name])['data']=posibility[x]
                s.append(graph.get_node(node_map[name]))                  
    temp=s
    s=[]
    for x in temp:
        graph.get_node(x['id'])['data']/=sumup
        s.append(graph.get_node(x['id']))
    ########################################################################
    #main program
    pagerank=PageRank(graph,s)
    p=pagerank.ApproximatePR()
    ########################################################################
    #form output
    for x in sorted(p,cmp=lambda x,y:cmp(p[y],p[x])):
        if p[x]>0 and reverse_node_map[x] in protein_name:
            print reverse_node_map[x],'\t',protein_name[reverse_node_map[x]],'\t',p[x]
        elif p[x]>0 and reverse_node_map[x] not in protein_name:
            print reverse_node_map[x],'\t','-','\t',p[x]
    # print graph.get_all_node_ids()
    # print len(s[0]['edges'])
