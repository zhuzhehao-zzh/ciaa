#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<random>
#include<cstdlib>
#include<ctime> 
#include<time.h> 
#include<chrono>
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;
using namespace std;

class Point
{
	private:
	    int pointId, clusterId;
	    
	    vector<double> values;

	public:
		int dimensions;
	    Point(int id,string line)
		{
	        dimensions=0;
	        pointId=id;
	        stringstream is(line);
	        double val;
	        while(is>>val)
			{
	            values.push_back(val);
	            dimensions++;
	        }
	        clusterId=0; 
	    }
	
	    int getDimensions()
		{
	        return dimensions;
	    }
	
	    int getCluster()
		{
	        return clusterId;
	    }
	
	    int getID()
		{
	        return pointId;
	    }
	
	    void setCluster(int val)
		{
	        clusterId=val;
	    }
	
	    double getVal(int pos)
		{
	        return values[pos];
	    }
};

class Cluster
{
	
	private:
	    int clusterId;
	    
	    vector<Point> points;
	
	public:
		vector<double> centroid;
		
	    Cluster(int clusterId,Point centroid)
		{
	        this->clusterId=clusterId;
	        for(int i=0;i<centroid.getDimensions();i++)
			{
	            this->centroid.push_back(centroid.getVal(i));
	        }
	        this->addPoint(centroid);
	    }
	
	    void addPoint(Point p)
		{
	        p.setCluster(this->clusterId);
	        points.push_back(p);
	    }
	
	    bool removePoint(int pointId)
		{
	        int size=points.size();
	
	        for(int i=0;i<size;i++)
	        {
	            if(points[i].getID()==pointId)
	            {
	                points.erase(points.begin()+i);
	                return true;
	            }
	        }
	        return false;
	    }
	
	    int getId()
		{
	        return clusterId;
	    }
	
	    Point getPoint(int pos)
		{
	        return points[pos];
	    }
	
	    int getSize()
		{
	        return points.size();
	    }
	
	    double getCentroidByPos(int pos)
	    {
	        return centroid[pos];
	    }
	
	    void setCentroidByPos(int pos,double val)
		{
	        this->centroid[pos]=val;
	    }
};

class CIAA{
	private:
	    int K, iters,dimensions,total_points;
	    vector<Cluster> clusters;

	    int getNearestClusterId(Point point,vector<double>& selfinf)
		{
	        double sum=0.0,min_dist;
	        int NearestClusterId;
	
	        for(int i=0;i<dimensions;i++)
	        {
	            sum+=pow(clusters[0].getCentroidByPos(i)-point.getVal(i), 2.0);
	        }
	
	        min_dist=sqrt(sum)+0*abs(selfinf[clusters[0].getId()]-selfinf[point.getID()]);
	        NearestClusterId=clusters[0].getId();
	
	        for(int i=1;i<K;i++)
	        {
	            double dist;
	            sum=0.0;
	
	            for(int j=0;j<dimensions;j++)
	            {
	                sum+=pow(clusters[i].getCentroidByPos(j)-point.getVal(j),2.0);
	            }
	
	            dist=sqrt(sum)+0*abs(selfinf[clusters[i].getId()]-selfinf[point.getID()]);
	
	            if(dist<min_dist)
	            {
	                min_dist=dist;
	                NearestClusterId=clusters[i].getId();
	            }
	        }
	
	        return NearestClusterId;
	    }
	
	public:
	    CIAA(int K,int iterations)
		{
	        this->K=K;
	        this->iters=iterations;
	    }
	
	    void run(vector<Point>& all_points,vector<double>& selfinf,vector<vector<int>>& edge)
		{
	        total_points=all_points.size();
	        dimensions=all_points[0].getDimensions();
	

	        vector<int> used_pointIds;
			//Initialize all the clusters
	        for(int i=1; i<=K; i++)
	        {
	            while(true)
	            {
	                int index=rand() % total_points;
	
	                if(find(used_pointIds.begin(),used_pointIds.end(), index)==used_pointIds.end())
	                {
	                    used_pointIds.push_back(index);
	                    all_points[index].setCluster(i);
	                    Cluster cluster(i,all_points[index]);
	                    clusters.push_back(cluster);
	                    break;
	                }
	            }
	        }
	        cout<<"Initial Clusters Size: "<<clusters.size()<<endl<<endl;
	
	
	        cout<<"CIAA Go"<<endl;
	
	        int iter=1;
	        while(true)
	        {
	            cout<<"Iteration Round is: "<<iter<<endl;
	            bool done=true;
	
	            for(int i=0;i<total_points;i++)
	            {
	                int currentClusterId=all_points[i].getCluster();
	                int nearestClusterId=getNearestClusterId(all_points[i],selfinf);
	
	                if(currentClusterId!=nearestClusterId)
	                {
	                    if(currentClusterId!=0)
						{
	                        for(int j=0;j<K;j++)
							{
	                            if(clusters[j].getId()==currentClusterId)
								{
	                                clusters[j].removePoint(all_points[i].getID());
	                            }
	                        }
	                    }
	
	                    for(int j=0;j<K;j++)
						{
	                        if(clusters[j].getId()==nearestClusterId)
							{
	                            clusters[j].addPoint(all_points[i]);
	                        }
	                    }
	                    all_points[i].setCluster(nearestClusterId);
	                    done = false;
	                }
	            }

	            for(int i=0;i<K;i++)
	            {
	                int ClusterSize=clusters[i].getSize();
					int min=999999,minmark=0;
					int cur=0;
					if (ClusterSize>0)
					{
						while (true)
						{
							if (cur>=ClusterSize) break;
							double sum=0,suminf=0;
							for (int t=0;t<ClusterSize;t++)
							{
								if (cur!=t)
								for (int j=0;j<dimensions;j++)
								{
									sum+=pow(clusters[i].getPoint(t).getVal(j)-clusters[i].getPoint(cur).getVal(j),2.0);
								}
								suminf+=0*abs(selfinf[clusters[i].getPoint(t).getID()]-selfinf[clusters[i].getPoint(cur).getID()]);
							}
							if (sqrt(sum)+suminf<min)
							{
								min=sqrt(sum);
								minmark=cur;
							}
							//cout<<clusters[i].getPoint(cur).getVal(0)<<" "<<clusters[i].getPoint(cur).getVal(1)<<" "<<sum<<endl;
 							cur++;
						}
					}
					//Update the centroids
					for (int j=0;j<dimensions;j++)
					{
						clusters[i].setCentroidByPos(j,clusters[i].getPoint(minmark).getVal(j));
					}
//					Check the Update
//	                for(int j = 0; j < dimensions; j++)
//	                {
//	                    double sum = 0.0;
//	                    if(ClusterSize > 0)
//	                    {
//	                        for(int p = 0; p < ClusterSize; p++)
//	                            sum += clusters[i].getPoint(p).getVal(j);
//	                        clusters[i].setCentroidByPos(j, sum / ClusterSize);
//	                    }
//	                }
	            }
	
	            if(done||iter>=iters)
	            {
	                cout<<"Clustering has finished in the round of: "<<iter<<endl<<endl;
	                break;
	            }
	            iter++;
	        }
	        //Check the Clusters
//	        int ans;
//			for (int i=0;i<K;i++)
//			{
//				for (int j=0;j<clusters[i].getSize();j++)
//				{
//					for (int k=0;k<dimensions;k++)
//					{				
//						cout<<clusters[i].getPoint(j).getVal(k)<<" ";
//					}
//					cout<<endl;
//					for (int k=0;k<dimensions;k++)
//					{
//						cout<<clusters[i].centroid[k]<<" ";
//						aver+=(clusters[i].getPoint(j).getVal(k)-clusters[i].centroid[k])*(clusters[i].getPoint(j).getVal(k)-clusters[i].centroid[k]);
//					}
//					cout<<endl;
//					for (int k=0;k<dimensions;k++)
//					{
//						cout<<(clusters[i].getPoint(j).getVal(k)-clusters[i].centroid[k])*(clusters[i].getPoint(j).getVal(k)-clusters[i].centroid[k])<<" ";
//						ans=ans+(clusters[i].getPoint(j).getVal(k)-clusters[i].centroid[k])*(clusters[i].getPoint(j).getVal(k)-clusters[i].centroid[k]);
//					}
//					cout<<endl;				 	 	 
//				}
//			}
//			cout<<ans<<endl; 

			double ans,sum;
			double aver,aver1;
			for (int i=0;i<K;i++)
			{
				for (int t=0;t<dimensions;t++)
				{
					aver=0;
					 
					for (int j=0;j<clusters[i].getSize();j++)
					{
						aver+=clusters[i].getPoint(j).getVal(t);
						aver1+=selfinf[clusters[i].getPoint(j).getID()];
					}
					aver=aver/clusters[i].getSize();
					aver1=aver1/clusters[i].getSize();
					sum=0;
					for (int j=0;j<clusters[i].getSize();j++)
					{
						sum+=pow(clusters[i].getPoint(j).getVal(t)-aver+0*(selfinf[clusters[i].getPoint(j).getID()]-aver1),2.0); 
					}
				}
				ans+=sum;
			}
			
			
			double sum1=0,sum2=0;
	        double sumd=0;
	        for(int i=0;i<K;i++)
	            {
	            	int ClusterSize=clusters[i].getSize();
	            	for (int t=0;t<ClusterSize;t++)
	            	{
	            		sumd=0;
	            		for (int j=0;j<dimensions;j++)
	            		{
	            			sumd+= pow(clusters[i].getCentroidByPos(j) - clusters[i].getPoint(t).getVal(j), 2.0);
						}
						sum1+=sqrt(sumd)+0*abs(selfinf[clusters[i].getPoint(t).getID()]-selfinf[clusters[i].getId()]);
					}
				}
				
			for(int i=0;i<K;i++)
	        {
	        	for (int t=0;t<K;t++)
	            {
	            	sumd=0;
	            	if (i!=t)
					for (int j=0;j<dimensions;j++)
	            	{
	            		sumd+= pow(clusters[i].getCentroidByPos(j) - clusters[t].getCentroidByPos(j), 2.0);
					}
					sum2+=sqrt(sumd)+0*abs(selfinf[clusters[t].getId()]-selfinf[clusters[i].getId()]);
				}
			}
			
			cout<<"Sparsity is: "<<(sum2/(K*(K-1)))/(sum1/all_points.size())<<endl;
			cout<<endl;
			
	        // Introduce the contents of different clusters
	        for(int i=0;i<K;i++){
	            cout<<"Cluster "<<clusters[i].getId()<<" : ";
	            for(int j=0;j<clusters[i].getSize();j++){
	                cout<<clusters[i].getPoint(j).getID()<<" ";
	            }
	            cout<<endl<<endl;
	        }
	        
	        int maxcluster=0, mincluster=9999;
	        double sumcluster;
	        for (int i=0;i<K;i++)
	        {
	        	if (clusters[i].getSize()>maxcluster) maxcluster=clusters[i].getSize();
	        	if (clusters[i].getSize()<mincluster) mincluster=clusters[i].getSize();
	        	sumcluster+=clusters[i].getSize();
			}
			sumcluster=sumcluster/K;
			cout<<"Cluster Distribution: "<<"MaxSize: "<<maxcluster<<" MinSize: "<<mincluster<<" AverageSize:"<<sumcluster<<endl;
			
			int maxinfgap=0,mininfgap=9999;
			double countinf=0,suminf=0;
			vector<double> infcom;
			for (int i=0;i<K;i++)
	        {
	        	if (clusters[i].getSize()==2)
	        	{
	                	if (abs(selfinf[clusters[i].getPoint(0).getID()]-clusters[i].getPoint(1).getID())>maxinfgap) 
							maxinfgap=abs(selfinf[clusters[i].getPoint(0).getID()]-clusters[i].getPoint(1).getID());
						if (abs(selfinf[clusters[i].getPoint(0).getID()]-clusters[i].getPoint(1).getID())<mininfgap) 
							mininfgap=abs(selfinf[clusters[i].getPoint(0).getID()]-clusters[i].getPoint(1).getID());
						infcom.push_back(abs(selfinf[clusters[i].getPoint(0).getID()]-clusters[i].getPoint(1).getID()));
						countinf++;
						suminf+=abs(selfinf[clusters[i].getPoint(0).getID()]-clusters[i].getPoint(1).getID());
	            }
				
			}
		
			//Dive inside the Clusters
			vector<vector<int>> checkk(all_points.size(),vector<int>(all_points.size(),0)); 
			int tier1=0,tier2=0,tier3=0;
			for (int i=0;i<K;i++)
			{
				for(int j=0;j<clusters[i].getSize();j++)
				{
					for(int t=0;t<clusters[i].getSize();t++)
					{
						if (j!=t&&abs(selfinf[clusters[i].getPoint(j).getID()]-clusters[i].getPoint(t).getID())<135)
						{
							if (checkk[clusters[i].getPoint(j).getID()][clusters[i].getPoint(t).getID()]==0)
							tier1++;
						}
						if (j!=t&&abs(selfinf[clusters[i].getPoint(j).getID()]-clusters[i].getPoint(t).getID())<270)
						{
							if (checkk[clusters[i].getPoint(j).getID()][clusters[i].getPoint(t).getID()]==0)
							tier2++;
						}
						if (j!=t&&abs(selfinf[clusters[i].getPoint(j).getID()]-clusters[i].getPoint(t).getID())<540)
						{
							if (checkk[clusters[i].getPoint(j).getID()][clusters[i].getPoint(t).getID()]==0)
							tier3++;
						}
						checkk[clusters[i].getPoint(j).getID()][clusters[i].getPoint(t).getID()]=1;
						checkk[clusters[i].getPoint(t).getID()][clusters[i].getPoint(j).getID()]=1;
						
					}
				}
			}
			cout<<"Num of different tiers: "<<tier1<<" "<<tier2<<" "<<tier3<<" "<<endl;
			
			
			sort(infcom.begin(),infcom.end());
			cout<<"Compare influence:"<<endl;
			//Especially check the results of clustes with only size 2 
			for (int i=0;i<infcom.size();i++)
				cout<<infcom[i]/2700<<" ";
			cout<<endl;
			
	        //Output cluster results
	        ofstream outfile;
	        outfile.open("clusters.txt");
	        if(outfile.is_open())
			{
	            for(int i=0; i<K; i++)
				{
	                for(int j=0; j<dimensions; j++)
					{
	                    cout<<clusters[i].getCentroidByPos(j)<<" ";
	                    outfile<<clusters[i].getCentroidByPos(j)<<" ";
	                }
	                cout<<endl;
	                outfile<<endl;
	            }
	            outfile.close();
	        }
	        else
			{
	            cout<<"Error: Unable to write to clusters.txt";
	        }
	
	    }
};

int interest[2700][10][5];
double interweight[10][5];
int NodeSize=2700,AreaSize=10,ActivitySize=5;

void BuildGraph(int size,vector<vector<int>>& edge,vector<int>& indegree,vector<int>& outdegree)
{
	unsigned seed;
    seed=time(0);
    srand(seed);
    
	int maxedge=size*(size-1)/2,minedge=size;
	int edgesize=rand()%(maxedge-minedge)+minedge;
	edgesize=5429;
	vector<int> f(size,0);
	int x,y,check;
	int edgecount=0;
	
	while (edgecount<edgesize)
	{
		x=rand()%size;
		if (f[x]==0)
		{
			
			y=rand()%size;
			if (x!=y && f[y]==0) 
			{
				//Check pair 
				//cout<<x<<" "<<y<<endl;
				edge[x][y]=1;
				f[x]=1;
				f[y]=1;
				outdegree[x]++;
				indegree[y]++;
				edgecount++;
			}
		}
		check=0;
		for (int i=0;i<size;i++)
			if (f[i]==0) check++;
		//cout<<check<<endl;
		if (check<2) break;
	}
	//cout<<edgesize<<" "<<edgecount<<endl;
	
	seed=time(0);
    srand(seed);
	int nx,ny;
	while (edgecount<edgesize)
	{
		x=rand()%size;
		y=rand()%size;
		//cout<<x<<" "<<y<<endl;
		if (x!=y&&edge[x][y]==0)
		{
			edge[x][y]=1;
			outdegree[x]++;
			indegree[y]++;
			edgecount++;
			if (edgecount==edgesize) break;
		}
	}
}

void quickSort(int left, int right, vector<vector<int>>& arr)
{
	if(left>=right) return;
	
	int i,j,base0,base1,base2,temp;
	i=left,j=right;
	
	base1=arr[left][1];
	base0=arr[left][0]; 
	base2=arr[left][2];
	
	while (i<j)
	{
		while (arr[j][1]>=base1 && i<j)
			j--;
		while (arr[i][1]<=base1 && i<j)
			i++;
		if(i<j)
		{
			temp=arr[i][1];
			arr[i][1]=arr[j][1];
			arr[j][1]=temp;
			
			temp=arr[i][0];
			arr[i][0]=arr[j][0];
			arr[j][0]=temp;
			
			temp=arr[i][2];
			arr[i][2]=arr[j][2];
			arr[j][2]=temp;
		}
	}
	
	arr[left][1]=arr[i][1];
	arr[i][1]=base1;
	
	arr[left][0]=arr[i][0];
	arr[i][0]=base0;
	
	arr[left][2]=arr[i][2];
	arr[i][2]=base2;
	
	
	quickSort(left,i-1,arr);
	quickSort(i+1,right,arr);
}

double Calmatch(int x, int y)
{
	double res;
	for (int i=0;i<AreaSize;i++)
	{
		for (int j=0;j<ActivitySize;j++)
		{
			if (interest[x][i][j]==interest[y][i][j])
				res+=interweight[i][j];
		} 
	}
	return res;
}

double Calmismatch(int x, int y)
{
	double res;
	for (int i=0;i<AreaSize;i++)
	{
		for (int j=0;j<ActivitySize;j++)
		{
			if (interest[x][i][j]!=interest[y][i][j])
				res+=interweight[i][j];
		} 
	}
	return res;
}

void Selfdiffusion (vector<vector<int>>& edge,vector<vector<int>>& node,vector<double>& selfinf)
{
	vector<double> tmpselfinf=selfinf;
	sort(tmpselfinf.begin(),tmpselfinf.end(),greater<int>());
	int threshold=tmpselfinf[NodeSize/10];	
} 

void CalSelf (vector<vector<int>>& edge,vector<vector<int>>& node,vector<double>& selfinf)
{
	int i,j;
	int size=node.size();
	for (i=0;i<size;i++)
		selfinf[i]=1;
		
//	Check the intial selfinf		
//	cout<<"The initial selfinf information:"<<endl;
//	for (i=0;i<size;i++)
//		cout<<selfinf[i]<<endl;

//	Check the in and out degrees		
//	for (i=0;i<size;i++)
//		cout<<node[i][0]<<" "<<node[i][1]<<" "<<node[i][2]<<endl;

	for (i=0;i<size;i++)
	{
		for (j=0;j<size;j++)
		{
			if (i!=j)
			{
				if (edge[node[i][0]][j]!=0)
				{
					selfinf[j]=selfinf[j]+(selfinf[node[i][0]])/node[i][2];
				}
			}
		}
	}
	cout<<"The selfinf information:"<<endl;
	cout<<size<<endl;
	int max=0;
	for (i=0;i<size;i++)
	{
		if (selfinf[i]>max) max=selfinf[i];
	}
	cout<<"maxinf:"<<max<<endl;
	
//  ReCheck Again
//	for (i=0;i<size;i++)
//		cout<<selfinf[i]<<endl;	
}
int main(){

    string filename="dataset_sample.txt";
    ifstream infile(filename.c_str());
	int K=10;
    int pointId=1;
    vector<Point> all_points;
    string line;

    while(getline(infile, line)){
        Point point(pointId, line);
        all_points.push_back(point);
        pointId++;
    }
    infile.close();
    cout<<"Successfully Read!"<<endl<<endl;
    if(all_points.size()<K){
        cout<<"Size error!"<<endl;
        return 1;
    }

    //CIAA go (Clustering with Influence Analysis)
    int iters = 50;
    
    vector<vector<int>> edge(all_points.size(),vector<int>(all_points.size())); 
    vector<int> indegree(all_points.size(),0);
    vector<int> outdegree(all_points.size(),0);
    
	vector<vector<int>> node(all_points.size(),vector<int>(3));
	vector<double> selfinf(all_points.size(),0);
	vector<vector<int>> interf(10,vector<int>(5));
	
	unsigned seed;
    seed=time(0);
    srand(seed);
    
    int tmpm,tmpn,tmpt;
	
	int intersum=0;
	for (int i=0;i<NodeSize;i++)
	{
		
		tmpt=rand()%AreaSize*ActivitySize;
		intersum+=tmpt;
		while (tmpt>0)
		{
			tmpm=rand()%AreaSize;
			tmpn=rand()%ActivitySize;
			if (interf[tmpm][tmpn]==0)
			{
				interest[i][tmpm][tmpn]=1;
				tmpt--;
			} 
		}
		
	}	
	//Calculate Weight
	for (int i=0;i<NodeSize;i++)
	{
		for (int j=0;j<AreaSize;j++)
		{
			for (int k=0;k<ActivitySize;k++)
			{
				interweight[j][k]+=interest[i][j][k];
			 } 
		}
	}
	for (int j=0;j<AreaSize;j++)
		{
			for (int k=0;k<ActivitySize;k++)
			{
				interweight[j][k]/=intersum;
			 } 
		}
    BuildGraph(all_points.size(),edge,indegree,outdegree);
    
//  Check the status of edges
//  	for (int i=0;i<all_points.size();i++)
//  		for (int j=0;j<all_points.size();j++)
//  			if (edge[i][j]==1) cout<<i<<" "<<j<<endl;
  			
  	for (int i=0;i<all_points.size();i++)
  	{
  		node[i][0]=i;
  		node[i][1]=indegree[i];
  		node[i][2]=outdegree[i];
	}
	quickSort(0,all_points.size()-1,node);
	
//	Check the node memory
//	cout<<"Check Size"<<node.size()<<" "<<all_points.size()<<endl;

// 	Check the sort result 
//	for (int i=0;i<all_points.size();i++)
//  	{
//		cout<<i<<" "<<indegree[i]<<" "<<outdegree[i]<<endl; 
//	}
//	cout<<"Here is sort results";
//	for (int i=0;i<all_points.size();i++)
//  	{
//		cout<<node[i][0]<<" "<<node[i][1]<<" "<<node[i][2]<<endl; 
//	}
//	for (int i=0;i<all_points.size();i++)
//  	{
//		cout<<nodeout[i][0]<<" "; 
//	}
//	cout<<endl;
//	for (int i=0;i<all_points.size();i++)
//  	{
//		cout<<nodeout[i][1]<<" "; 
//	}
//	cout<<endl;
//	for (int i=0;i<all_points.size();i++)
//  	{
//		cout<<nodeout[i][2]<<" "; 
//	}
//	cout<<endl;
	
	CalSelf(edge,node,selfinf);
	
	int sumd,sumd_t=0,sumi=0,max=0;
	for (int i=0;i<all_points.size();i++)
		for (int j=0;j<all_points.size();j++)
		if (i!=j)
		{
			sumd_t=0;
			for (int t=0;t<all_points[0].dimensions;t++)
				sumd_t+=pow(all_points[i].getVal(t)-all_points[j].getVal(t),2.0);
			if (sumd_t>max) max=sumd_t; 
			sumd+=sqrt(sumd_t);
		}
	cout<<"maxdist:"<<max<<endl; 
	
    high_resolution_clock::time_point beginTime = high_resolution_clock::now();
   	CIAA ciaa(K,iters);
	ciaa.run(all_points,selfinf,edge);
	
	high_resolution_clock::time_point endTime = high_resolution_clock::now();
    microseconds timeInterval1 = std::chrono::duration_cast<microseconds>(endTime - beginTime);
    
	cout<<"²åÈëºÄÊ±£º"<<timeInterval1.count()<<"Î¢Ãë"<<endl;
	cout<<"CIAA done"<<endl;
	
    return 0;
}
