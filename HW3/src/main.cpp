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
	float x, y;
};

struct block{
	int id;
	int area;
	bool rotate;
	int w, h;
	float x, y;
};

struct net{
	int deg;
	vector<int> blockList, terminalList;
};

struct node{
	int id;
	int idx;
	int parent = -5;
	int lc = -5, rc = -5;
	int w = 0, h = 0;
	int x, y;
	vector<tuple<int, int, int, int> > whpair;
};

//global var
int numBlock, numTerminal, numNet, numPin;
int total_area, boundary, wire_length;
int cur_w, cur_h;
unordered_map<int, block*> blocks;
unordered_map<int, terminal*> terminals;
vector<net> nets;
node* initTree, *bestTree, *localTree;
int m1,m2,m3;
//

void read_block(FILE* input){
	
	fscanf(input, "NumHardRectilinearBlocks : %d\n", &numBlock);
	fscanf(input, "NumTerminals : %d\n\n", &numTerminal);
	int cnt = numBlock;
	while(cnt--){
	
		int a;
		block* tmp = new block();
		fscanf(input, " sb%d hardrectilinear 4 (%f, %f) (%d, %d) (%d, %d) (%d, %d)\n", &tmp->id, &tmp->x, &tmp->y, &a, &a, &tmp->w, &tmp->h, &a, &a);
	
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
		fscanf(input, " p%d %f %f\n", &tmp->id, &tmp->x, &tmp->y);
		terminals[tmp->id] = tmp;
	}
	return;
}
void read_net(FILE* input){
	fscanf(input, " NumNets : %d\n", &numNet);
	fscanf(input, " NumPins : %d\n", &numPin);
	int cnt = numNet;
	while(cnt--){
		net tmp;
		string str;
		char c[30];
		fscanf(input, " NetDegree : %d\n", &tmp.deg);		
		for(int i=0; i<tmp.deg; i++){
			fscanf(input, " %s\n", c);
			str = c;
			if(str[0] == 'p'){
				str.erase(0,1);
				tmp.terminalList.push_back(stoi(str));
			}else{
				str.erase(0,2);
				tmp.blockList.push_back(stoi(str));
			}
		}
		nets.push_back(tmp);
	}
	
	return;
}

vector<int> init_PE(){
	//form PE = 12V3V4V...nV
	vector<int> PE;
	PE.push_back(0);
	for(int i=1;i<numBlock;i++){
		PE.push_back(i);
		int m = rand()%2;
		if(m)PE.push_back(-1);
		else PE.push_back(-1);
	}
	return PE;
}

//build initTree from PE

void build_tree(vector<int> PE){
	stack<int> stk;
	initTree = new node[2*numBlock-1]();
	int idx1, idx2;
	//use idx to indicate parent and children
	for(int i=0; i<PE.size(); i++){
		initTree[i].id = PE[i];
		if(PE[i]>=0){
			stk.push(i);
			initTree[i].w = blocks[PE[i]]->w;
			initTree[i].h = blocks[PE[i]]->h;
			//cout<<PE[i]<<" "<<initTree[i].w<<" "<<initTree[i].h<<" ";
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
		}else{
			cout<<"error build tree\n";
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
		rc_x = tree[idx].x + tree[tree[idx].lc].w;
		rc_y = tree[idx].y;
		cal_xy(tree, tree[idx].lc, lc_x, lc_y);
		cal_xy(tree, tree[idx].rc, rc_x, rc_y);
	}else if(tree[idx].id == -2){//H
		lc_x = tree[idx].x;
		lc_y = tree[idx].y;
		rc_x = tree[idx].x;
		rc_y = tree[idx].y + tree[tree[idx].lc].h;
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
		//cout<<sboth.size()<<" " << tmp << "   ";
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
		}else if(tree[tmp].id == -2){//H
		
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
		}else{
			cout<<tree[tmp].id<<" "<<tree[tree[tmp].parent].id<<" "<<tree[tmp].lc<<" "<<tree[tmp].rc<<endl;
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
//	cost += A;
	if(cur_h > boundary){
		cost += (cur_h - boundary)*(cur_h - boundary);
	}
	if(cur_w > boundary){
		cost += (cur_w - boundary)*(cur_w - boundary);
	}
//	for(auto it = blocks.begin();it!= blocks.end(); it++){
//		if(it->second->x > boundary){
//			cost += 1000;
//		}
//		if(it->second->y > boundary){
//			cost += 1000;
//		}
//	}
	//wire length
	wire_length = 0;
	for(int i=0; i<nets.size();i++){
		float maxx = 0, maxy = 0;
		float minx = 1000000, miny = 1000000;
		float tmpx, tmpy;
		for(auto j=nets[i].blockList.begin(); j!=nets[i].blockList.end(); j++){
			tmpx = blocks[*j]->x + blocks[*j]->w/2.0;
			tmpy = blocks[*j]->y + blocks[*j]->h/2.0;
			maxx = tmpx >maxx? tmpx:maxx;
			maxy = tmpy >maxy? tmpy:maxy;
			minx = tmpx <minx? tmpx:minx;
			miny = tmpy <miny? tmpy:miny;
		}
		for(auto j=nets[i].terminalList.begin(); j!=nets[i].terminalList.end(); j++){
			tmpx = terminals[*j]->x;
			tmpy = terminals[*j]->y;
			maxx = tmpx >maxx? tmpx:maxx;
			maxy = tmpy >maxy? tmpy:maxy;
			minx = tmpx <minx? tmpx:minx;
			miny = tmpy <miny? tmpy:miny;
		}
		wire_length += (maxx-minx)+(maxy-miny);
	}
	cost += (int)(0.01 * wire_length);
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
	if(idx1 == -1 || idx2 == -1) cout<<"idx == -1"<<endl;
	if(PE[idx1]<0 ||PE[idx2]<0) cout<<"PE[idx]<0"<<endl;
	//change order in PE
	int tmp = PE[idx1];
	PE[idx1]= PE[idx2];
	PE[idx2] = tmp; 
	//change parent
	tmp = tree[idx1].parent;
	tree[idx1].parent = tree[idx2].parent;
	tree[idx2].parent = tmp;
	//change order in tree
	node n;
	n = tree[idx1];
	tree[idx1] = tree[idx2];
	tree[idx2] = n;
	
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
	while(!s1.empty() && !s2.empty()){
		if(s1.top() == s2.top()){
			sboth.push(s1.top());
			s1.pop();
			s2.pop();
		}else{
			break;
		}
	}
	while(!s1.empty()){
		sboth.push(s1.top());
		s1.pop();
	}
	while(!s2.empty()){
		sboth.push(s2.top());
		s2.pop();
	}
	area_computation(sboth, tree, PE);

	return PE;
} 

vector<int> M2_complement(vector<int>PE, node* tree){
	int cnt=0, idx, tmp;
	int picked = rand() % (numBlock - 1) + 1;//range(0, number of 'V' and 'H')
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
	while(PE[idx-1]<0){
		idx = idx-1;
	}
	while(PE[idx] == -1 || PE[idx] == -2){ //HVH or VHV, progress once if it's not a chain
		tree[idx].id = tree[idx].id==-1? -2:-1;
		PE[idx] = PE[idx]==-1? -2:-1;
		stk.push(idx);
		if(idx == PE.size()-1) break;
		if(PE[idx+1] == -1 ||PE[idx+1] == -2) idx = idx+1;
		else break;
	}
	while(tree[idx].parent != -5 &&(tree[tree[idx].parent].id == -1 || tree[tree[idx].parent].id == -2)){
		stk.push(tree[idx].parent);
		idx = tree[idx].parent;
		//cout<<idx<<" ";
	}

//	idx = tmp;
//	
//	while(tree[idx].id == -1 || tree[idx].id == -2){
//		stk.push(idx);
//		idx = tree[idx].parent;
//		if(idx == -5) break;
//	}
	while(!stk.empty()){
		sboth.push(stk.top());
		stk.pop();
	}

	area_computation(sboth, tree, PE);
	return PE;
}

vector<int> M3_swap(vector<int>PE, node* tree){
	int picked, min=0, opcnt, cnt, idx1, idx2, terminate=0, tmp;
	bool done = false;
	while(!done){
		if(terminate++>100){
			return PE;
		}
		opcnt = cnt = 0;
		picked = rand() % (numBlock) + 1; //1~numBlock
		for(int i=0; i<PE.size(); i++){
			if(PE[i]==-1 || PE[i]==-2){
				opcnt++;
			}else if(PE[i]>=0){
				cnt++;
				if(cnt == picked){
					//cout<<"0.5";
//					if(PE[i+1]<0 && (PE[i+1] != PE[i-1]) && (cnt-1>opcnt+1)){//change operand i and operator i+1
//						idx1 = i;
//						idx2 = i+1;
//						done = true;
//						
//						tmp = PE[idx1];
//						PE[idx1]= PE[idx2];
//						PE[idx2] = tmp;
//						cout<<"1";
//						//change tree
//						tmp = tree[tree[idx2].parent].lc;
//						tree[tree[idx2].parent].lc = idx1;
//						tree[tree[idx2].parent].rc = idx2;
//						tree[idx1].parent = tree[idx2].parent;
//						tree[idx2].rc = tree[idx2].lc;
//						tree[idx2].lc = tmp;
//						tree[tmp].parent = idx1;
//						tree[tree[idx2].rc].parent = idx1;
//						
//						cout<<"3";
//						node n;
//						n = tree[idx1];
//						tree[idx1] = tree[idx2];
//						tree[idx2] = n;
//						
//						
//					}else 
					if(i-1>0 && PE[i-1]<0 && (PE[i+1] != PE[i-1])){//change operand i and operator i-1
						idx1 = i;
						idx2 = i-1;
						done = true;
						//change order in PE
						tmp = PE[idx1];
						PE[idx1]= PE[idx2];
						PE[idx2] = tmp;
						//change H/V's left child's parent to its parent
						tmp = tree[idx2].lc;
						tree[tree[idx2].lc].parent = tree[idx2].parent;
						tree[tree[idx2].parent].lc = tree[idx2].lc;
						//change H/V's left child to right child
						tree[tree[idx2].rc].parent = idx1;
						tree[idx2].lc = tree[idx2].rc;	
						tree[idx2].parent = tree[idx1].parent;
						tree[idx2].rc = idx2;
						tree[idx1].parent = idx1;
						//change order in tree
						node n;
						n = tree[idx1];
						tree[idx1] = tree[idx2];
						tree[idx2] = n;
					} 
					break;
				}
			}
		}
	}
	node n;
	//change the area vector on path
	stack<int> s1, s2, sboth;
	n = tree[idx1];
	while(n.parent!=-5){
		if(tree[n.parent].id >=0) cout<<"error 1 "<<tree[n.parent].id<<" "<<n.parent<<endl;
		s1.push(n.parent);
		n = tree[n.parent];
	}
	n = tree[idx2]; //idx1 is now H/V's idx
	while(n.parent!=-5){
		if(tree[n.parent].id >=0) cout<<"error 2 "<<tree[n.parent].id<<" "<<n.parent<<endl;
		s2.push(n.parent);
		n = tree[n.parent];
	}
	while(!s1.empty() && !s2.empty()){
		if(s1.top() == s2.top()){
			sboth.push(s1.top());
			s1.pop();
			s2.pop();
		}else{
			break;
		}
	}
	while(!s1.empty()){
		sboth.push(s1.top());
		s1.pop();
	}
	while(!s2.empty()){
		sboth.push(s2.top());
		s2.pop();
	}
	area_computation(sboth, tree, PE);	

	m3++;
	return PE;
}

void copy_tree(node* t1, node* t2, int k){	
	for(int i=0; i<k; i++){
		t1[i] = t2[i];
	} 
	return;
}

int cal_initT(float p, vector<int> PE, node* tree){
	node* tmp_tree;
	tmp_tree = new node[2*numBlock-1]();
	vector<int> tmp_PE = PE;
	copy_tree(tmp_tree, tree, PE.size());
	
	int tmp_cost = cal_cost();
	long long move, avg=0, tmp;
	int cnt = 1000;
	while(cnt-- >0){
		move = rand() %3;
		//tmp_PE = M2_complement(tmp_PE, t);
		if(move == 0){
			tmp_PE = M1_swap(tmp_PE, tmp_tree); 
		}else if(move == 1){
			tmp_PE = M2_complement(tmp_PE, tmp_tree); 
		}else if(move == 2){
			tmp_PE = M3_swap(tmp_PE, tmp_tree); 
		}
		tmp = cal_cost();
		if(tmp > tmp_cost){
			avg += (tmp - tmp_cost);
		}
		tmp_cost = tmp;
	}
	cout<<(double)(1.0*avg/1000)<<endl;
	return (double)(1.0*avg/1000)/(-1*log(p));
	
}

vector<int> SA_floorplanning(vector<int> PE, node* tree, float r, int k, float p, int e){// init sol, tree, rate for cooling, var k, accepting rate p
	vector<int> E = PE, NE;
	vector<int> Best = E;
	stack<int> emptystk;
	bestTree = new node[2*numBlock-1]();
	localTree = new node[2*numBlock-1]();
	
	copy_tree(bestTree, tree, PE.size());
	copy_tree(localTree, tree, PE.size());
	area_computation(emptystk, bestTree, E);
	int T0 = cal_initT(p, PE, tree);
	cout<<"init T0: "<<T0<<endl;
	double T =T0; 
	
	int M = 0, MT = 0, uphill = 0, reject =0;
	int N = k*numBlock;
	
	copy_tree(tree, bestTree, PE.size());
	area_computation(emptystk, bestTree, E);
	int cost = cal_cost(), d_cost = 0, n_cost, b_cost = cost;
	
	while(true){
		MT = uphill = reject = 0;
		while(true){
			int move = rand() %3;
			if(move == 0){
				cout<<"1";
				m1++;
				NE = M1_swap(E, tree); 
				cout<<"\r";
			}else if(move == 1){
				cout<<"2";
				m2++;
				NE = M2_complement(E, tree); 
				cout<<"\r";
			}else if(move == 2){
				cout<<"3";
				NE = M3_swap(E, tree); 
				cout<<"\r";
			}
			MT++;
			n_cost = cal_cost();//cost of NE
			d_cost =  n_cost - cost;//cost of E
			cout<<"T: "<<T<<" w: "<<cur_w<<" h: "<<cur_h<<" wire length: "<<wire_length<<" cost: "<<n_cost<<"\r";
			//cout<<"done move"<<endl;
			//cout<<"exp: "<<exp(-1*d_cost/T)<<"\r";
			if(d_cost<=0 || (double) rand() / (RAND_MAX + 1.0) < exp(-1.0*d_cost/T)){
				//cout<<"accept ans:"<<d_cost<<endl;
				if(d_cost>0) uphill++;
				E = NE;
				copy_tree(localTree, tree, PE.size());
				cost = n_cost;
				if(cost<b_cost){
					Best = E;
					b_cost = cost;
					copy_tree(bestTree, tree, PE.size());
				}	
			}else{
					reject++;
					copy_tree(tree, localTree, PE.size());
			}
//			cout<<endl;
//			int aaaaaa;
//			for(int i = 0; i<199;i++){
//				if(E[i] >= 0){
//					cout<<"("<<Best[i]<<")";
//				}else if(Best[i] == -1){
//					cout<<"V";
//				}else if(Best[i] == -2){
//					cout<<"H";
//				}
//			}
//			cout << endl;
//			cin>>aaaaaa;
			if(uphill>N || MT>2*N) break;
		}
		T = r*T;
//		area_computation(emptystk, bestTree, E);
//		for(int i = 0; i<199;i++){
//			if(Best[i] >= 0){
//				cout<<"("<<Best[i]<<")";
//			}else if(Best[i] == -1){
//				cout<<"V";
//			}else if(Best[i] == -2){
//				cout<<"H";
//			}
//		}
		cout<<"T: "<<T<<" w: "<<cur_w<<" h: "<<cur_h<<" wire length: "<<wire_length<<" " << (1.0*reject/MT) <<" "<< endl;
		if((1.0*reject/MT) > 0.95 || T<e){
			break;	
		} 
		
	}
	
	return Best;
} 

int main(void){
	FILE* input_block = fopen("../testcase/n100.hardblocks", "r");
	FILE* input_net = fopen("../testcase/n100.nets", "r");
	FILE* input_pl = fopen("../testcase/n100.pl", "r");
	srand(time(NULL));
	float r_deadspace = 0.15;
	read_block(input_block);
	read_net(input_net);
	read_terminal(input_pl);
	vector<int> PE = init_PE();
	build_tree(PE);
	stack<int> emptystk;
	area_computation(emptystk, initTree, PE);
	boundary = (int)sqrt(total_area*(1+r_deadspace));
	cout<<"Boundary: "<<boundary<<endl;
//	bestTree = new node[2*numBlock-1]();
//	copy_tree(bestTree, initTree, 199);
	PE = SA_floorplanning(PE, initTree, 0.85, 5, 0.9, 1000);
//	node * cmp_tree = new node[2*numBlock-1]();
//	for(int i = 0;i<50;i++){
//		PE = M2_complement(PE,initTree);
//	}
//
//	for(int i = 0;i<100;i++){
//		
//		PE = M3_swap(PE,initTree);
//	}
//	copy_tree(cmp_tree, initTree, PE.size());
//	build_tree(PE);
//	cout<<"===============parent============\n";
//	
//	for(int i = 0;i<PE.size();i++){
//		if(cmp_tree[i].parent != initTree[i].parent){
//			cout<<"error "<<i<<" "<<cmp_tree[i].parent<<" "<<initTree[i].parent;
//		}
//	}
//	cout<<"\n=================lc==============\n";
//	
//	for(int i = 0;i<PE.size();i++){
//		if(cmp_tree[i].lc!= initTree[i].lc){
//			cout<<"error "<<i<<" "<<cmp_tree[i].lc<<" "<<initTree[i].lc;
//		}
//	}
//	cout<<"\n=================rc==============\n";
//	
//	for(int i = 0;i<PE.size();i++){
//		if(cmp_tree[i].rc != initTree[i].rc){
//			cout<<"error "<<i<<" "<<cmp_tree[i].rc<<" "<<initTree[i].rc;
//		}
//	}
//	cout<<"\n================id===============\n";
//	for(int i = 0;i<PE.size();i++){
//		if(cmp_tree[i].id != initTree[i].id){
//			cout<<"error "<<i<<" "<<cmp_tree[i].id<<" "<<initTree[i].id;
//		}
//	}
//	cout<<"=================================\n";
//	
//	area_computation(emptystk, bestTree, PE);
//	if(cur_w < boundary && cur_h < boundary){
//		cout<<"success"<<endl;
//	}
//	cout<<cur_w<<" "<<cur_h<<" "<<wire_length;
	//build_tree(PE);
	//bestTree = new node[199]();
	//copy_tree(bestTree, initTree, 199); 
	int opcnt=0, cnt=0;
	for(int i=0;i<PE.size();i++){
		if(initTree[i].id == -1){
			cout<<"V";	
			opcnt++;
		}
		else if(initTree[i].id == -2){
			cout<<"H";
			opcnt++;	
		} 
		else{
			cout<<"("<<initTree[i].id<<")";	
			cnt++;
		} 
		if(opcnt>cnt) cout<<"illegal\n";
	}

	cout<<endl<<m1<<" "<<m2<<" "<<m3<<endl;
	cout<<total_area;
	
	//area_computation(emptystk, initTree, PE);
	ofstream output;
	output.open("./out.txt");
	output<<"Wirelength "<<wire_length<<endl;
	output<<"Blocks\n";
	for(int i=0;i<numBlock;i++){
		//cout<<"sb"<<i<<" "<<blocks[i]->w<<" "<<blocks[i]->h<<" "<<blocks[i]->rotate<<endl;
		output<<"sb"<<blocks[i]->id<<" "<<blocks[i]->x<<" "<<blocks[i]->y<<" "<<blocks[i]->rotate<<endl;
	} 
	output.close();
	return 0;
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
//	for(int i=0;i<PE.size();i++){
//		if(initTree[i].id == -1) cout<<"V";
//		else if(initTree[i].id == -2) cout<<"H";
//		else cout<<"("<<initTree[i].id<<")";
//	}
//	cout<<endl;
//	
