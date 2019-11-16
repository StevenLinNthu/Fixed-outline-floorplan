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

//build tree from PE
node* tree;
void build_tree(vector<int> PE){
	stack<int> stk;
	tree = new node[2*numBlock]();
	int idx1, idx2;
	//use idx to indicate parent and children
	for(int i=0; i<PE.size(); i++){
		tree[i].id = PE[i];
		if(PE[i]>=0){
			stk.push(i);
			tree[i].w = blocks[PE[i]]->w;
			tree[i].h = blocks[PE[i]]->h;
			tree[i].whpair.push_back(make_tuple(tree[i].w, tree[i].h, -1, -1));
			if(tree[i].w != tree[i].h){
				tree[i].whpair.push_back(make_tuple(tree[i].h, tree[i].w,-1, -1));
			} 
		}else if(PE[i]==-1){ //'V'
			idx2 = tree[i].rc = stk.top();
			stk.pop();
			idx1 = tree[i].lc = stk.top();
			stk.pop();
			stk.push(i);
			tree[tree[i].rc].parent = tree[tree[i].lc].parent = i;
			for(int j=0; j<tree[idx1].whpair.size(); j++){
				for(int k=0; k<tree[idx2].whpair.size(); k++){
					int w=0, h=0;
					w = get<0>(tree[idx1].whpair[j]) + get<0>(tree[idx2].whpair[k]);
					h = max(get<1>(tree[idx1].whpair[j]), get<1>(tree[idx2].whpair[k]));
					if(tree[i].w==0 && tree[i].h==0){
						tree[i].w = w;
						tree[i].h = h;
						tree[i].whpair.push_back(make_tuple(w, h, j, k));
					}else if(w < tree[i].w || h< tree[i].h){
						tree[i].whpair.push_back(make_tuple(w, h, j, k));
						tree[i].w = tree[i].w>w? w:tree[i].w;
						tree[i].h = tree[i].h>h? h:tree[i].h;
					}
				}
			}
		}else if(PE[i]==-2){ //'H'
			idx2 = tree[i].rc = stk.top();
			stk.pop();
			idx1 = tree[i].lc = stk.top();
			stk.pop();
			stk.push(i);
			tree[tree[i].rc].parent = tree[tree[i].lc].parent = i;
			for(int j=0; j<tree[idx1].whpair.size(); j++){
				for(int k=0; k<tree[idx2].whpair.size(); k++){
					int w=0, h=0;
					w = max(get<0>(tree[idx1].whpair[j]), get<0>(tree[idx2].whpair[k]));
					h = get<1>(tree[idx1].whpair[j]) + get<1>(tree[idx2].whpair[k]);
					if(tree[i].w==0 && tree[i].h==0){
						tree[i].w = w;
						tree[i].h = h;
						tree[i].whpair.push_back(make_tuple(w, h, j, k));
					}else if(w < tree[i].w || h< tree[i].h){
						tree[i].whpair.push_back(make_tuple(w, h, j, k));
						tree[i].w = tree[i].w>w? w:tree[i].w;
						tree[i].h = tree[i].h>h? h:tree[i].h;
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
//	int min = 1000000;
//	for(int i = 0; i<tree[PE.size()-1].whpair.size();i++){
//		int idx, w, h;
//		w = get<0>(tree[PE.size()-1].whpair[i]);
//		h = get<1>(tree[PE.size()-1].whpair[i]);
//		if(w*h<min){
//			min = w*h;
//			cout<<w<<" "<<h<<endl;
//		}
//	}
	
	return 0;
}
