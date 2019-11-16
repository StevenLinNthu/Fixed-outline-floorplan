#include<iostream>
#include<unistd.h>
#include<cstdlib>
#include<fstream>
#include<cstdio>
#include<string>
#include<stack>
#include<vector>
#include<tuple>
#include<list>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<ctime>
#include<cstdlib>
#include<cstring>
#include<sys/time.h>
#include<unordered_map>

using namespace std;


struct terminal{
	int id;
	int x, y;
};

struct block{
	int id;
	int area;
	bool rotate;
	int w, h, x, y;
	int parent, lc, rc;
	//pointer?
};

struct net{
	int deg;
	vector<int> blockList, terminalList;
};

struct node{
	int id;
	int idx;
	int parent=-5, lc=-5, rc=-5;//parent, left child, right child
	int w=0, h=0;
	vector<tuple<int, int, int, int> > whpair;
};

//global var
int numBlock, numTerminal, numNet, numPin;
int total_area;
unordered_map<int, block*> blocks;
unordered_map<int, terminal*> terminals;
unordered_map<int, net*> nets;
node* initTree, bestTree, localTree;
//

void read_block(FILE* input){
	
	fscanf(input, "NumHardRectilinearBlocks : %d\n", &numBlock);
	fscanf(input, "NumTerminals : %d\n\n", &numTerminal);
	int cnt = numBlock;
	while(cnt--){
	
		int a;
		block* tmp = new block();
		fscanf(input, " sb%d hardrectilinear 4 (%d, %d) (%d, %d) (%d, %d) (%d, %d)\n", &tmp->id, &tmp->x, &tmp->y, &a, &a, &tmp->w, &tmp->h, &a, &a);
	
		tmp->area = tmp->w * tmp->h;
		total_area += tmp->area;
		tmp->rotate = false;
		blocks[tmp->id] = tmp;
	}
	return;
}

void read_terminal(FILE* input){
	int cnt = numTerminal;
	while(cnt--){
		terminal* tmp = new terminal();
		fscanf(input, " p%d %d %d\n", &tmp->id, &tmp->x, &tmp->y);
		terminals[tmp->id] = tmp;
	}
	return;
}
void read_net(FILE* input){
	fscanf(input, " NumNets : %d\n", &numNet);
	fscanf(input, " NumPins : %d\n", &numPin);
	int cnt = numNet;
	while(cnt--){
		net* tmp = new net();
		string str;
		char c[30];
		fscanf(input, " NetDegree : %d\n", &tmp->deg);		
		for(int i=0; i<tmp->deg; i++){
			fscanf(input, " %s\n", c);
			str = c;
			if(str[0] == 'p'){
				str.erase(0,1);
				tmp->terminalList.push_back(stoi(str));
			}else{
				str.erase(0,2);
				tmp->blockList.push_back(stoi(str));
			}
		}
	}
	return;
}

vector<int> init_PE(){
	//form PE = 12V3V4V...nV
	vector<int> PE;
	unordered_map<int,block*>::iterator it=blocks.begin();
	PE.push_back(it->second->id);
	it++;
	PE.push_back(it->second->id);
	it++;
	while(it!=blocks.end()){
		PE.push_back(-1);
		PE.push_back(it->second->id);
		it++;
	}
	PE.push_back(-1);	
	return PE;
}

//build initTree from PE

void build_tree(vector<int> PE){
	stack<int> stk;
	initTree = new node[2*numBlock]();
	int idx1, idx2;
	//use idx to indicate parent and children
	for(int i=0; i<PE.size(); i++){
		initTree[i].id = PE[i];
		if(PE[i]>=0){
			stk.push(i);
			initTree[i].w = blocks[PE[i]]->w;
			initTree[i].h = blocks[PE[i]]->h;
			initTree[i].whpair.push_back(make_tuple(initTree[i].w, initTree[i].h, -1, -1));
			if(initTree[i].w != initTree[i].h){
				initTree[i].whpair.push_back(make_tuple(initTree[i].h, initTree[i].w,-1, -1));
			} 
		}else if(PE[i]==-1){ //'V'
			idx2 = initTree[i].rc = stk.top();
			stk.pop();
			idx1 = initTree[i].lc = stk.top();
			stk.pop();
			stk.push(i);
			initTree[initTree[i].rc].parent = initTree[initTree[i].lc].parent = i;
			for(int j=0; j<initTree[idx1].whpair.size(); j++){
				for(int k=0; k<initTree[idx2].whpair.size(); k++){
					int w=0, h=0;
					w = get<0>(initTree[idx1].whpair[j]) + get<0>(initTree[idx2].whpair[k]);
					h = max(get<1>(initTree[idx1].whpair[j]), get<1>(initTree[idx2].whpair[k]));
					if(initTree[i].w==0 && initTree[i].h==0){
						initTree[i].w = w;
						initTree[i].h = h;
						initTree[i].whpair.push_back(make_tuple(w, h, j, k));
					}else if(w < initTree[i].w || h< initTree[i].h){
						initTree[i].whpair.push_back(make_tuple(w, h, j, k));
						initTree[i].w = initTree[i].w>w? w:initTree[i].w;
						initTree[i].h = initTree[i].h>h? h:initTree[i].h;
					}
				}
			}
		}else if(PE[i]==-2){ //'H'
			idx2 = initTree[i].rc = stk.top();
			stk.pop();
			idx1 = initTree[i].lc = stk.top();
			stk.pop();
			stk.push(i);
			initTree[initTree[i].rc].parent = initTree[initTree[i].lc].parent = i;
			for(int j=0; j<initTree[idx1].whpair.size(); j++){
				for(int k=0; k<initTree[idx2].whpair.size(); k++){
					int w=0, h=0;
					w = max(get<0>(initTree[idx1].whpair[j]), get<0>(initTree[idx2].whpair[k]));
					h = get<1>(initTree[idx1].whpair[j]) + get<1>(initTree[idx2].whpair[k]);
					if(initTree[i].w==0 && initTree[i].h==0){
						initTree[i].w = w;
						initTree[i].h = h;
						initTree[i].whpair.push_back(make_tuple(w, h, j, k));
					}else if(w < initTree[i].w || h< initTree[i].h){
						initTree[i].whpair.push_back(make_tuple(w, h, j, k));
						initTree[i].w = initTree[i].w>w? w:initTree[i].w;
						initTree[i].h = initTree[i].h>h? h:initTree[i].h;
					}
				}
			}
		}
	}
	return;
}

//M1 swap two adjacent operands
void M1_swap(vector<int> PE, node* tree){
	//randomly choose a number i between 0 and numBlock-1, swap the ith operands in PE with i+1th operands
	int min=0, cnt=0;
	int picked = rand() % (numBlock - min + 1) + min;
	int idx1 = -1, idx2 = -1;	
	for(int i=0; i<PE.size(); i++){
		if(PE[i]>=0){
			cnt++;
			if(cnt == picked){
				idx1 = i;
			}
			if(idx1 >= 0){
				idx2 = i;
				break;
			}
		}
	}
	//change order in PE
	int tmp = PE[idx1];
	PE[idx1]= PE[idx2];
	PE[idx2] = tmp; 
	//change order in tree
	node n;
	n = tree[idx1];
	tree[idx1] = tree[idx2];
	tree[idx2] = n;
	//change parent
	tmp = tree[idx1].parent;
	tree[idx1].parent = tree[idx2].parent;
	tree[idx2].parent = tmp;
	//change the area vector on path
	stack<int> s1, s2, sboth;
	n = tree[idx1];
	while(n.parent!=-5){
		s1.push(n.parent);
		n = tree[n.parent];
	}
	n = tree[idx2];
	while(n.parent!=-5){
		s2.push(n.parent);
		n = tree[n.parent];
	}
	while(true){
		if(s1.top() == s2.top()){
			sboth.push(s1.top());
			s1.pop();
			s2.pop();
		}else if(!s1.empty()){
			sboth.push(s1.top());
			s1.pop();
		}else if(!s2.empty()){
			sboth.push(s2.top());
			s2.pop();
		}else{
			break;
		}
	}
	while(!sboth.empty()){
		tmp = sboth.top();
		sboth.pop();
		if(tree[tmp].id == -1){//V
			idx2 = tree[tmp].rc;
			idx1 = tree[tmp].lc;
			tree[tmp].w = tree[tmp].h = 0;
			tree[tmp].whpair.clear(); 
			for(int j=0; j<tree[idx1].whpair.size(); j++){
				for(int k=0; k<tree[idx2].whpair.size(); k++){
					int w=0, h=0;
					w = get<0>(tree[idx1].whpair[j]) + get<0>(tree[idx2].whpair[k]);
					h = max(get<1>(tree[idx1].whpair[j]), get<1>(tree[idx2].whpair[k]));
					if(tree[tmp].w==0 && tree[tmp].h==0){
						tree[tmp].w = w;
						tree[tmp].h = h;
						tree[tmp].whpair.push_back(make_tuple(w, h, j, k));
					}else if(w < tree[tmp].w || h< tree[tmp].h){
						tree[tmp].whpair.push_back(make_tuple(w, h, j, k));
						tree[tmp].w = tree[tmp].w>w? w:tree[tmp].w;
						tree[tmp].h = tree[tmp].h>h? h:tree[tmp].h;
					}
				}
			}
		}else{//H
			idx2 = tree[tmp].rc;
			idx1 = tree[tmp].lc;
			tree[tmp].w = tree[tmp].h = 0;
			tree[tmp].whpair.clear(); 
			for(int j=0; j<tree[idx1].whpair.size(); j++){
				for(int k=0; k<tree[idx2].whpair.size(); k++){
					int w=0, h=0;
					w = max(get<0>(tree[idx1].whpair[j]), get<0>(tree[idx2].whpair[k]));
					h = get<1>(tree[idx1].whpair[j]) + get<1>(tree[idx2].whpair[k]);
					if(tree[tmp].w==0 && tree[tmp].h==0){
						tree[tmp].w = w;
						tree[tmp].h = h;
						tree[tmp].whpair.push_back(make_tuple(w, h, j, k));
					}else if(w < tree[tmp].w || h< tree[tmp].h){
						tree[tmp].whpair.push_back(make_tuple(w, h, j, k));
						tree[tmp].w = tree[tmp].w>w? w:tree[tmp].w;
						tree[tmp].h = tree[tmp].h>h? h:tree[tmp].h;
					}
				}
			}
		}
	}	
	return;
} 
int SA_floorplanning(int k){
//	E = init_sol();
//	Best = E;
//	float T0 = 10000;
//	int M = 0, MT = 0, uphill = 0;
//	int N = k*numBlock;
//	while(...){
//		
//	}
	
	
	return 1;
} 

int main(void){
	FILE* input_block = fopen("../testcase/n100.hardblocks", "r");
	read_block(input_block);
	vector<int> PE = init_PE();
	build_tree(PE);
//pick smallest area pair from root
//	int min = 1000000;
//	for(int i = 0; i<initTree[PE.size()-1].whpair.size();i++){
//		int idx, w, h;
//		w = get<0>(initTree[PE.size()-1].whpair[i]);
//		h = get<1>(initTree[PE.size()-1].whpair[i]);
//		if(w*h<min){
//			min = w*h;
//			cout<<w<<" "<<h<<endl;
//		}
//	}
	
	return 0;
}
