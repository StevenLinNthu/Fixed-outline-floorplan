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
	int x, y;
	vector<tuple<int, int, int, int> > whpair;
};

//global var
int numBlock, numTerminal, numNet, numPin;
int total_area, boundary, wire_length;
int cur_w, cur_h;
unordered_map<int, block*> blocks;
unordered_map<int, terminal*> terminals;
unordered_map<int, net*> nets;
node* initTree, *bestTree, *localTree;
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
void cal_wh(node* tree, int picked, int idx){
	int lc_idx, rc_idx, lc_picked, rc_picked;
	tree[idx].w = get<0>(tree[idx].whpair[picked]);
	tree[idx].h = get<1>(tree[idx].whpair[picked]);
	lc_picked = get<2>(tree[idx].whpair[picked]);
	rc_picked = get<3>(tree[idx].whpair[picked]);
	lc_idx = tree[idx].lc;
	rc_idx = tree[idx].rc;
	if(lc_idx != -5) cal_wh(tree, lc_picked, lc_idx);
	if(rc_idx != -5) cal_wh(tree, rc_picked, rc_idx);
	return;
}
void cal_xy(node* tree, int idx, int x, int y){
	tree[idx].x = x;
	tree[idx].y = y;
	int lc_x, lc_y, rc_x, rc_y;
	if(tree[idx].id == -1){//V
		lc_x = tree[idx].x;
		lc_y = tree[idx].y;
		rc_x = tree[idx].x + tree[tree[idx].lc].x;
		rc_y = tree[idx].y;
		cal_xy(tree, tree[idx].lc, lc_x, lc_y);
		cal_xy(tree, tree[idx].rc, rc_x, rc_y);
	}else if(tree[idx].id == -2){//H
		lc_x = tree[idx].x;
		lc_y = tree[idx].y;
		rc_x = tree[idx].x;
		rc_y = tree[idx].y + tree[tree[idx].lc].y;
		cal_xy(tree, tree[idx].lc, lc_x, lc_y);
		cal_xy(tree, tree[idx].rc, rc_x, rc_y);
	}else if(tree[idx].id >=0){
		if(blocks[tree[idx].id]->w != tree[idx].w){
			int tmp;
			tmp = blocks[tree[idx].id]->w;
			blocks[tree[idx].id]->w = blocks[tree[idx].id]->h;
			blocks[tree[idx].id]->h = tmp;
			blocks[tree[idx].id]->rotate = !blocks[tree[idx].id]->rotate;
		}
		blocks[tree[idx].id]->x = x;
		blocks[tree[idx].id]->y = y;
	}
	return;
}

void area_computation(stack<int> sboth, node* tree, vector<int> PE){
	int tmp, idx1, idx2;
	
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
	//compute x, y coordinate
	int min = 100000000, root_idx = PE.size()-1;
	int lc_idx, rc_idx, lc_picked, rc_picked;
	int idx, w, h, picked;
	for(int i = 0; i<tree[root_idx].whpair.size();i++){
		
		w = get<0>(tree[root_idx].whpair[i]);
		h = get<1>(tree[root_idx].whpair[i]);
		if(w*h<min){
			min = w*h;
			cur_w = w;
			cur_h = h;
			picked = i;
		}
	}
	cal_wh(tree, picked, root_idx);
	cal_xy(tree, root_idx, 0, 0);
	return;
}

int cal_cost(){
	int cost = 0;
	int A = cur_h * cur_w;
	cost += A;
	if(cur_h > boundary && cur_w > boundary){
		cost += (A - boundary*boundary)*50000;
	}else if(cur_h > boundary){
		cost += (cur_h - boundary)*cur_w*50000;
	}else if(cur_w > boundary){
		cost += (cur_w - boundary)*cur_h*50000;
	}
	//wire length
	wire_length = 0;
	for(auto i=nets.begin(); i!=nets.end();i++){
		int maxx = 0, maxy = 0;
		int minx = 1000000, miny = 1000000;
		int tmpx, tmpy;
		for(auto j=i->second->blockList.begin(); j!=i->second->blockList.end(); j++){
			tmpx = blocks[*j]->x + blocks[*j]->w/2;
			tmpy = blocks[*j]->y + blocks[*j]->h/2;
			maxx = tmpx >maxx? tmpx:maxx;
			maxy = tmpy >maxy? tmpy:maxy;
			maxx = tmpx <minx? tmpx:minx;
			maxy = tmpy <miny? tmpy:miny;
		}
		for(auto j=i->second->terminalList.begin(); j!=i->second->terminalList.end(); j++){
			tmpx = terminals[*j]->x;
			tmpy = terminals[*j]->y;
			maxx = tmpx >maxx? tmpx:maxx;
			maxy = tmpy >maxy? tmpy:maxy;
			maxx = tmpx <minx? tmpx:minx;
			maxy = tmpy <miny? tmpy:miny;
		}
		wire_length += (maxx-minx)+(maxy-miny);
	}
	cost += wire_length;
	return cost;
}

//M1 swap two adjacent operands
vector<int> M1_swap(vector<int> PE, node* tree){
	//randomly choose a number i between 0 and numBlock-1, swap the ith operands in PE with i+1th operands
	int min=1, cnt=0;
	int picked = rand() % (numBlock - min) + min, picked2 = picked;
	while(picked == picked2){
		picked2 = rand() % (numBlock - min) + min;
	}
	int idx1 = -1, idx2 = -1;	
	for(int i=0; i<PE.size(); i++){
		if(PE[i]>=0){
			cnt++;
			if(cnt == picked){
				idx1 = i;
			}
			if(cnt == picked2){
				idx2 = i;
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
	//cout<<"2"<<endl;
	//change parent
	tmp = tree[idx1].parent;
	tree[idx1].parent = tree[idx2].parent;
	tree[idx2].parent = tmp;
	//change the area vector on path
	//cout<<"3"<<endl;
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
	while(s1.empty() && s2.empty()){
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
		}
	}
	area_computation(sboth, tree, PE);	
	return PE;
} 

vector<int> M2_complement(vector<int>PE, node* tree){
	int min=1, cnt=0, idx, tmp;
	int picked = rand() % (PE.size()-numBlock - min) + min;//range(0, number of 'V' and 'H')
	stack<int> stk, sboth;
	for(int i=0; i<PE.size(); i++){
		if(PE[i]<0){
			cnt++;
			if(cnt == picked){
				idx = i;
				break;
			}
		}
	}
	
	while(tree[idx].rc != -5 && (tree[tree[idx].rc].id == -1 || tree[tree[idx].rc].id == -2)){//find the chain content VHVHVHV....
		picked = tree[idx].rc;
	}
	tmp = idx;
	while(tree[idx].id == -1 || tree[idx].id == -2){ //HVH or VHV, progress once if it's not a chain
		tree[idx].id = tree[idx].id==-1? -2:-1;
		if(tree[idx].parent == -5) break;
		if(idx == tree[tree[idx].parent].rc) idx = tree[idx].parent;
		else break;
	}
	idx = tmp;
	while(tree[idx].id == -1 || tree[idx].id == -2){
		stk.push(idx);
		idx = tree[idx].parent;
		if(idx == -5) break;
	}
	while(!stk.empty()){
		sboth.push(stk.top());
		stk.pop();
	}
	area_computation(sboth, tree, PE);	
	return PE;
}

vector<int> M3_swap(vector<int>PE, node* tree){
	int picked, min=0, opcnt, cnt, idx1, idx2, terminate=0;
	bool done = false;
	while(!done){
		if(terminate++<1000){
			return PE;
		}
		opcnt = cnt = 0;
		picked = rand() % (numBlock - min + 1) + min;
		for(int i=0; i<PE.size(); i++){
			if(PE[i]==-1 || PE[i]==-2){
				opcnt++;
			}else if(PE[i]>=0){
				cnt++;
				if(cnt == picked){
					//cout<<"0.5";
					if(PE[i+1]<0 && PE[i+1] != PE[i-1] && 2*(opcnt+1)<i){
						cout<<"0.6";
						idx1 = i;
						idx2 = i+1;
						done = true;
					} 
					break;
				}
			}
		}
	}
	cout<<"1"<<endl;
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
	while(s1.empty() && s2.empty()){
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
		}
	}
	area_computation(sboth, tree, PE);	
	return PE;
}

void copy_tree(node* t1, node* t2, int k){	
	for(int i=0; i<k; i++){
		t1[i] = t2[i];
	} 
}

int cal_initT(float p, vector<int> PE, node* tree){
	node tmpTree[PE.size()];
	node* t = tmpTree;
	vector<int> tmp_PE = PE;
	copy_tree(t, tree, PE.size());
	int tmp_cost = cal_cost();
	int move, avg, tmp;
	int cnt = 1000;
	while(cnt-- >0){
		move = rand() %3;
		//tmp_PE = M2_complement(tmp_PE, t);
		if(move == 0){
			tmp_PE = M1_swap(tmp_PE, t); 
		}else if(move == 1){
			tmp_PE = M2_complement(tmp_PE, t); 
		}else if(move == 2){
			tmp_PE = M3_swap(tmp_PE, t); 
		}
		tmp = cal_cost();
		avg += tmp - tmp_cost;
		tmp_cost = tmp;
	}
	return (int)(avg/1000)/log(p);
	
}

vector<int> SA_floorplanning(vector<int> PE, node* tree, float r, int k, float p, int e){// init sol, tree, rate for cooling, var k, accepting rate p
	vector<int> E = PE, NE;
	vector<int> Best = E;
	bestTree = new node[2*numBlock]();
	localTree = new node[2*numBlock]();
	copy_tree(bestTree, tree, PE.size());
	copy_tree(localTree, tree, PE.size());
	int T0 = cal_initT(p, PE, tree);
	double T =T0; 
	
	int M = 0, MT = 0, uphill = 0, reject =0;
	int N = k*numBlock;
	stack<int> emptystk;
	copy_tree(tree, bestTree, PE.size());
	area_computation(emptystk, bestTree, E);
	int cost = cal_cost(), d_cost = 0, n_cost, b_cost = cost;
	
	while(true){
		MT = uphill = reject = 0;
		while(true){
			int move = rand() %3;
			if(move == 0){
				//cout<<"1";
				NE = M1_swap(E, tree); 
			}else if(move == 1){
				//cout<<"2";
				NE = M2_complement(E, tree); 
			}else if(move == 2){
				//cout<<"3";
				NE = M3_swap(E, tree); 
			}
			MT++;
			n_cost = cal_cost();//cost of NE
			d_cost =  n_cost - cost;//cost of E
			//cout<<"done move"<<endl;
			if(d_cost<=0 || (double) rand() / (RAND_MAX + 1.0) < exp(-1*d_cost/T)){
				//cout<<"accept ans:"<<d_cost<<endl;
				if(d_cost>0) uphill++;
				E = NE;
				copy_tree(localTree, tree, PE.size());
				cost = n_cost;
				if(cost<b_cost){
					Best = E;
					copy_tree(tree, bestTree, PE.size());
				}	
			}else{
					reject++;
					copy_tree(tree, localTree, PE.size());
			}
			if(uphill>N || MT>2*N) break;
		}
		T = r*T;
		cout<<"T: "<<T<<endl;
		if((reject/MT) > 0.95 || T<e){
			break;	
		} 
		
	}
	cout<<"done";
	
	return Best;
} 

int main(void){
	FILE* input_block = fopen("../testcase/n100.hardblocks", "r");
	FILE* input_net = fopen("../testcase/n100.nets", "r");
	FILE* input_pl = fopen("../testcase/n100.pl", "r");
	float r_deadspace = 0.15;
	read_block(input_block);
	read_net(input_net);
	read_terminal(input_pl);
	vector<int> PE = init_PE();
	build_tree(PE);
	stack<int> emptystk;
	area_computation(emptystk, initTree, PE);
	
	boundary = (int)sqrt(total_area*(1+r_deadspace));
	PE = SA_floorplanning(PE, initTree, 0.9, 5, 0.9, 10);
	area_computation(emptystk, bestTree, PE);
	if(cur_w < boundary && cur_h < boundary){
		cout<<"success"<<endl;
	}
	//test
//	vector<int> PE;
//	PE.push_back(1);
//	PE.push_back(2);
//	PE.push_back(-2);
//	PE.push_back(3);
//	PE.push_back(4);
//	PE.push_back(-1);
//	PE.push_back(5);
//	PE.push_back(6);
//	PE.push_back(-1);
//	PE.push_back(-2);
//	PE.push_back(-1);
//	build_tree(PE);
//	PE = M3_swap(PE, initTree);
	for(int i=0;i<PE.size();i++){
		if(initTree[i].id == -1) cout<<"V";
		else if(initTree[i].id == -2) cout<<"H";
		else cout<<"("<<initTree[i].id<<")";
	}
	//M2_complement(PE, initTree);
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
/**
done(11)(98)V(97)V(96)V(45)V(44)V(43)V(42)V(41)V(40)V(39)V(38)V(37)V(36)V(35)V(34)V(33)V(32)V(31)V(30)V(29)V(28)V(27)V(26)V(25)V(24)V(83)V(22)V(77)V(8)V(7)V(6)V(5)V(4)V(3)V(2)V(1)V(0)V(10)V(99)V(12)V(13)V(14)V(15)V(16)V(17)V(18)V(19)V(20)V(21)V(46)V(47)V(48)V(49)V(50)V(51)V(9)V(53)V(54)V(55)V(56)V(57)V(58)V(59)V(60)V(61)V(62)V(63)V(64)V(65)V(66)V(67)V(68)V(69)V(70)V(71)V(72)V(73)V(74)V(75)V(76)V(52)V(78)H(79)V(80)V(81)V(82)V(23)V(84)V(85)V(86)V(87)V(88)V(89)V(90)V(91)V(92)V(93)V(94)V(95)V
**/
