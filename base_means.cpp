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
	    int dimensions;
	    vector<double> values;

	public:
	    Point(int id, string line)
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
	
	    void setCentroidByPos(int pos, double val)
		{
	        this->centroid[pos]=val;
	    }
};

class KMeans{
	private:
	    int K,iters,dimensions,total_points;
	    vector<Cluster> clusters;

	    int getNearestClusterId(Point point)
		{
	        double sum=0.0, min_dist;
	        int NearestClusterId;
	
	        for(int i=0;i<dimensions;i++)
	        {
	            sum+=pow(clusters[0].getCentroidByPos(i)-point.getVal(i), 2.0);
	        }
	
	        min_dist=sqrt(sum);
	        NearestClusterId=clusters[0].getId();
	
	        for(int i=1;i<K;i++)
	        {
	            double dist;
	            sum=0.0;
	
	            for(int j=0;j<dimensions;j++)
	            {
	                sum+=pow(clusters[i].getCentroidByPos(j)-point.getVal(j), 2.0);
	            }
	
	            dist=sqrt(sum);
	
	            if(dist<min_dist)
	            {
	                min_dist=dist;
	                NearestClusterId=clusters[i].getId();
	            }
	        }
	
	        return NearestClusterId;
	    }
	
	public:
	    KMeans(int K, int iterations)
		{
	        this->K=K;
	        this->iters=iterations;
	    }
	
	    void run(vector<Point>& all_points
		){
	
	        total_points=all_points.size();
	        dimensions=all_points[0].getDimensions();
	
	        vector<int> used_pointIds;
	
	        for(int i=1;i<=K;i++)
	        {
	            while(true)
	            {
	                int index=rand()%total_points;
	
	                if(find(used_pointIds.begin(),used_pointIds.end(),index)==used_pointIds.end())
	                {
	                    used_pointIds.push_back(index);
	                    all_points[index].setCluster(i);
	                    Cluster cluster(i, all_points[index]);
	                    clusters.push_back(cluster);
	                    break;
	                }
	            }
	        }
 			cout<<"Initial Clusters Size: "<<clusters.size()<<endl<<endl;
	
	
	        cout<<"K-means Go"<<endl;
	
	        int iter=1;
	        while(true)
	        {
	            cout<<"Iteration Round is: "<<iter<<endl;
	            bool done = true;
	            for(int i = 0; i < total_points; i++)
	            {
	                int currentClusterId=all_points[i].getCluster();
	                int nearestClusterId=getNearestClusterId(all_points[i]);
	
	                if(currentClusterId!=nearestClusterId)
	                {
	                    if(currentClusterId!=0)
						{
	                        for(int j=0; j<K; j++)
							{
	                            if(clusters[j].getId()==currentClusterId){
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
	                    done=false;
	                }
	            }
	            for(int i=0;i<K;i++)
	            {
	                int ClusterSize=clusters[i].getSize();
	
	                for(int j=0;j<dimensions;j++)
	                {
	                    double sum=0.0;
	                    if(ClusterSize>0)
	                    {
	                        for(int p=0;p<ClusterSize;p++)
	                            sum+=clusters[i].getPoint(p).getVal(j);
	                        clusters[i].setCentroidByPos(j,sum/ClusterSize);
	                    }
	                }
	            }
	
	            if(done||iter>=iters)
	            {
	                cout<<"Clustering has finished in the round of: "<<iter<<endl<<endl;
	                break;
	            }
	            iter++;
	        }
	        for(int i=0; i<K; i++)
			{
	            cout<<"Cluster "<<clusters[i].getId()<<" : ";
	            for(int j=0; j<clusters[i].getSize(); j++){
	                cout<<clusters[i].getPoint(j).getID()<<" ";
	            }
	            cout<<endl<<endl;
	        }
	        
			double ans,sum;
			double aver;
			for (int i=0;i<K;i++)
			{
				for (int t=0;t<dimensions;t++)
				{
					aver=0;
					for (int j=0;j<clusters[i].getSize();j++)
					{
						aver+=clusters[i].getPoint(j).getVal(t);
					}
					aver=aver/clusters[i].getSize();
					
					sum=0;
					for (int j=0;j<clusters[i].getSize();j++)
					{
						sum+=pow(clusters[i].getPoint(j).getVal(t)-aver,2.0); 
					}
				}
				ans+=sum;
			}
			cout<<ans<<endl;
		}
};
void BuildGraph(int size,vector<vector<int>>& edge,vector<int>& edgenum)
{
	unsigned seed;
    seed=time(0);
    srand(seed);
    
	int maxedge=size*(size-1)/2,minedge=size;
	int edgesize=rand()%(maxedge-minedge)+minedge;
	
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
				edge[x][y]=1;
				f[x]=1;
				f[y]=1;
				edgenum[x]++;
				edgenum[y]++;
				edgecount++;
			}
		}
		check=0;
		for (int i=0;i<size;i++)
			if (f[i]==0) check++;
		if (check<2) break;
	}
	seed=time(0);
    srand(seed);
	while (edgecount<edgesize)
	{
		x=rand()%size;
		y=rand()%size;
		if (x!=y&&edge[x][y]==0)
		{
			edge[x][y]=1;
			edgenum[x]++;
			edgenum[y]++;
			edgecount++;
			if (edgecount==edgesize) break;
		}
	}
}
void quickSort(int left, int right, vector<vector<int>>& arr)
{
	if(left>=right) return;
	
	int i,j,base,baset,temp;
	i=left,j=right;
	base=arr[left][0];
	baset=arr[left][1]; 
	while (i<j)
	{
		while (arr[j][0]>=base && i<j)
			j--;
		while (arr[i][0]<=base && i<j)
			i++;
		if(i<j)
		{
			temp=arr[i][0];
			arr[i][0]=arr[j][0];
			arr[j][0]=temp;
			
			temp=arr[i][1];
			arr[i][1]=arr[j][1];
			arr[j][1]=temp;
		}
	}
	arr[left][0]=arr[i][0];
	arr[i][0]=base;
	
	arr[left][1]=arr[i][1];
	arr[i][1]=baset;
	
	quickSort(left,i-1,arr);
	quickSort(i+1,right,arr);
}
int main()
{
    string filename="dataset_Iris.txt";
    ifstream infile(filename.c_str());
	int K=5;
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
    int iters=500;
    vector<vector<int>> edge(all_points.size(),vector<int>(all_points.size())); 
    vector<int> edgenum(all_points.size(),0);
	vector<vector<int>> nodeout(all_points.size(),vector<int>(2));
	vector<int> selfinf(all_points.size(),0);
	
    BuildGraph(all_points.size(),edge,edgenum);
//  	for (int i=0;i<all_points.size();i++)
//  		for (int j=0;j<all_points.size();j++)
//  			if (edge[i][j]==1) cout<<i<<" "<<j<<endl;
  	for (int i=0;i<all_points.size();i++)
  	{
  		nodeout[i][0]=edgenum[i];
  		nodeout[i][1]=i;
	}
	quickSort(0,all_points.size()-1,nodeout);
	high_resolution_clock::time_point beginTime = high_resolution_clock::now();
	KMeans kmeans(K,iters);
    kmeans.run(all_points);
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    microseconds timeInterval1 = std::chrono::duration_cast<microseconds>(endTime - beginTime);
	cout<<"²åÈëºÄÊ±£º"<<timeInterval1.count()<<"Î¢Ãë"<<endl;
    return 0;
}
