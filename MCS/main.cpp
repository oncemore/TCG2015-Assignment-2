/*
 * This code is provied as a sample code of Hw 2 of "Theory of Computer Game".
 * The "genmove" function will randomly output one of the legal move.
 * This code can only be used within the class.
 *
 * 2015 Nov. Hung-Jui Chang
 * */
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>

#define BOARDSIZE        9
#define BOUNDARYSIZE    11
#define COMMANDLENGTH 1000
#define DEFAULTTIME     10
#define DEFAULTKOMI      7

#define SAMPLES         50

#define MAXGAMELENGTH 1000
#define MAXSTRING       50
#define MAXDIRECTION     4

#define NUMINTERSECTION 81
#define HISTORYLENGTH   200

#define EMPTY            0
#define BLACK            1
#define WHITE            2
#define BOUNDARY         3

#define SELF             1
#define OPPONENT         2

#define NUMGTPCOMMANDS      15

#define LOCALVERSION      1
#define GTPVERSION        2

#define RANDOM           1
#define NONRANDOM        0

using namespace std;
int _board_size = BOARDSIZE;
int _board_boundary = BOUNDARYSIZE;
double _komi =  DEFAULTKOMI;
const int DirectionX[MAXDIRECTION] = {-1, 0, 1, 0};
const int DirectionY[MAXDIRECTION] = { 0, 1, 0,-1};
const char LabelX[]="0ABCDEFGHJ";

void do_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int move);
double final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]);

// MCST tree
struct node {
  struct node *firstchild;
  struct node *nextsibling;
  int branch_id;
  int games;
  float score;
};
/*
 * This function reset the board, the board intersections are labeled with 0,
 * the boundary intersections are labeled with 3.
 * */
void reset(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    for (int i = 1 ; i <= BOARDSIZE; ++i) {
	for (int j = 1 ; j <= BOARDSIZE; ++j) {
	    Board[i][j] = EMPTY;
	}
    }
    for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
	Board[0][i] = Board[BOUNDARYSIZE-1][i] = Board[i][0] = Board[i][BOUNDARYSIZE-1] = BOUNDARY;
    }
}

/*
 * This function return the total number of liberity of the string of (X, Y) and
 * the string will be label with 'label'.
 * */
int find_liberty(int X, int Y, int label, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE]) {
    // Label the current intersection
    ConnectBoard[X][Y] |= label;
    int total_liberty = 0;
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
	// Check this intersection has been visited or not
	if ((ConnectBoard[X+DirectionX[d]][Y+DirectionY[d]] & (1<<label) )!= 0) continue;

	// Check this intersection is not visited yet
	ConnectBoard[X+DirectionX[d]][Y+DirectionY[d]] |=(1<<label);
	// This neighboorhood is empty
	if (Board[X+DirectionX[d]][Y+DirectionY[d]] == EMPTY){
	    total_liberty++;
	}
	if(total_liberty >= 2)
        break;
	// This neighboorhood is self stone
	else if (Board[X+DirectionX[d]][Y+DirectionY[d]] == Board[X][Y]) {
	    total_liberty += find_liberty(X+DirectionX[d], Y+DirectionY[d], label, Board, ConnectBoard);
	}
	// We only consider the case that liberty is 1 or 0.
	if(total_liberty >= 2)
        break;
    }
    return total_liberty;
}

/*
 * This function count the liberties of the given intersection's neighboorhod
 * */
void count_liberty(int X, int Y, int Board[BOUNDARYSIZE][BOUNDARYSIZE], int Liberties[MAXDIRECTION]) {
    int ConnectBoard[BOUNDARYSIZE][BOUNDARYSIZE];
    // Initial the ConnectBoard
    for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
        for (int j = 0 ; j < BOUNDARYSIZE; ++j) {
            ConnectBoard[i][j] = 0;
        }
    }
    // Find the same connect component and its liberity
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
        Liberties[d] = 0;
        if (Board[X+DirectionX[d]][Y+DirectionY[d]] == BLACK ||
            Board[X+DirectionX[d]][Y+DirectionY[d]] == WHITE    ) {
            Liberties[d] = find_liberty(X+DirectionX[d], Y+DirectionY[d], d, Board, ConnectBoard);
        }
    }
}

/*
 * This function count the number of empty, self, opponent, and boundary intersections of the neighboorhood
 * and saves the type in NeighboorhoodState.
 * */
void count_neighboorhood_state(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn, int* empt, int* self, int* oppo ,int* boun, int NeighboorhoodState[MAXDIRECTION]) {
    for (int d = 0 ; d < MAXDIRECTION; ++d) {
	// check the number of nonempty neighbor
	switch(Board[X+DirectionX[d]][Y+DirectionY[d]]) {
	    case EMPTY:    (*empt)++;
			   NeighboorhoodState[d] = EMPTY;
			   break;
	    case BLACK:    if (turn == BLACK) {
			       (*self)++;
			       NeighboorhoodState[d] = SELF;
			   }
			   else {
			       (*oppo)++;
			       NeighboorhoodState[d] = OPPONENT;
			   }
			   break;
	    case WHITE:    if (turn == WHITE) {
			       (*self)++;
			       NeighboorhoodState[d] = SELF;
			   }
			   else {
			       (*oppo)++;
			       NeighboorhoodState[d] = OPPONENT;
			   }
			   break;
	    case BOUNDARY: (*boun)++;
			   NeighboorhoodState[d] = BOUNDARY;
			   break;
	}
    }
}

/*
 * This function remove the connect component contains (X, Y) with color "turn"
 * And return the number of remove stones.
 * */
int remove_piece(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    int remove_stones = (Board[X][Y]==EMPTY)?0:1;
    Board[X][Y] = EMPTY;
    for (int d = 0; d < MAXDIRECTION; ++d) {
        if (Board[X+DirectionX[d]][Y+DirectionY[d]] == turn) {
            remove_stones += remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], turn);
        }
    }
    return remove_stones;
}
/*
 * This function update Board with place turn's piece at (X,Y).
 * Note that this function will not check if (X, Y) is a legal move or not.
 * */
void update_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int Liberties[4];
    int NeighboorhoodState[4];
    count_neighboorhood_state(Board, X, Y, turn,
	    &num_neighborhood_empt,
	    &num_neighborhood_self,
	    &num_neighborhood_oppo,
	    &num_neighborhood_boun, NeighboorhoodState);
    // check if there is opponent piece in the neighboorhood
    if (num_neighborhood_oppo != 0) {

    count_liberty(X, Y, Board, Liberties);

	for (int d = 0 ; d < MAXDIRECTION; ++d) {
	    // check if there is opponent component only one liberty
	    if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1 && Board[X+DirectionX[d]][Y+DirectionY[d]]!=EMPTY) {
            remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], Board[X+DirectionX[d]][Y+DirectionY[d]]);
	    }
	}
    }
    Board[X][Y] = turn;
}
/*
 * This function update Board with place turn's piece at (X,Y).
 * Note that this function will check if (X, Y) is a legal move or not.
 * */
int update_board_check(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int X, int Y, int turn) {
    // Check the given coordination is legal or not
    if ( X < 1 || X > BOARDSIZE || Y < 1 || Y > BOARDSIZE || Board[X][Y]!=EMPTY)
	return 0;
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int Liberties[4];
    int NeighboorhoodState[4];
    count_neighboorhood_state(Board, X, Y, turn,
	    &num_neighborhood_empt,
	    &num_neighborhood_self,
	    &num_neighborhood_oppo,
	    &num_neighborhood_boun, NeighboorhoodState);
    // Check if the move is a legal move
    // Condition 1: there is a empty intersection in the neighboorhood
    int legal_flag = 0;
    count_liberty(X, Y, Board, Liberties);
    if (num_neighborhood_empt != 0) {
	legal_flag = 1;
    }
    else {
	// Condition 2: there is a self string has more than one liberty
	for (int d = 0; d < MAXDIRECTION; ++d) {
	    if (NeighboorhoodState[d] == SELF && Liberties[d] > 1) {
		legal_flag = 1;
	    }
	}
	if (legal_flag == 0) {
	// Condition 3: there is a opponent string has exactly one liberty
	    for (int d = 0; d < MAXDIRECTION; ++d) {
		if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1) {
		    legal_flag = 1;
		}
	    }
	}
    }

    if (legal_flag == 1) {
    // check if there is opponent piece in the neighboorhood
	if (num_neighborhood_oppo != 0) {
	    for (int d = 0 ; d < MAXDIRECTION; ++d) {
		// check if there is opponent component only one liberty
		if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1 && Board[X+DirectionX[d]][Y+DirectionY[d]]!=EMPTY) {
		    remove_piece(Board, X+DirectionX[d], Y+DirectionY[d], Board[X+DirectionX[d]][Y+DirectionY[d]]);
		}
	    }
	}
	Board[X][Y] = turn;
    }

    return (legal_flag==1)?1:0;
}

/*
 * This function return the number of legal moves with color "turn" and
 * saves all legal moves in MoveList
 * */
int gen_legal_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length
                   , int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE], int MoveList[HISTORYLENGTH], int random_bit) {
    int NextBoard[BOUNDARYSIZE][BOUNDARYSIZE];
    int num_neighborhood_self = 0;
    int num_neighborhood_oppo = 0;
    int num_neighborhood_empt = 0;
    int num_neighborhood_boun = 0;
    int legal_moves = 0;
    int x,y;
    int next_x, next_y;
    int random_index;
    int Liberties[4];
    int liberty_counted = 0;
    int NeighboorhoodState[4];
    bool eat_move = 0;

    int intersection_checked = 0;

    random_index = rand()%81;

    while(intersection_checked<81){
        if(random_bit == 1){
            x = random_index / 9 + 1;
            y = random_index % 9 + 1;
            random_index = (random_index + 67) % 81;
        }
        else{
            x = (intersection_checked / 9) + 1;
            y = (intersection_checked % 9) + 1;
        }
	    // check if current
	    if (Board[x][y] == EMPTY) {
		// check the liberty of the neighborhood intersections
		num_neighborhood_self = 0;
		num_neighborhood_oppo = 0;
		num_neighborhood_empt = 0;
		num_neighborhood_boun = 0;
		// count the number of empty, self, opponent, and boundary neighboorhood
		count_neighboorhood_state(Board, x, y, turn,
			&num_neighborhood_empt,
			&num_neighborhood_self,
			&num_neighborhood_oppo,
			&num_neighborhood_boun, NeighboorhoodState);
        liberty_counted = 0;
		// check if the empty intersection is a legal move
		if(num_neighborhood_empt == 4){
            next_x = x;
            next_y = y;
            eat_move = 0;
		}
		else {
            next_x = next_y = 0;
            eat_move = 0;
            count_liberty(x, y, Board, Liberties);
            // Case 1: exist empty intersection in the neighborhood
             if (num_neighborhood_empt > 0) {
                 next_x = x;
                 next_y = y;
                 // check if it is a capture move
                 for (int d = 0 ; d < MAXDIRECTION; ++d) {
                     if (NeighboorhoodState[d] == OPPONENT && Liberties[d] == 1) {
                         eat_move = 1;
                     }
                 }
             }
             // Case 2: no empty intersection in the neighborhood
             else {
                // Case 2.1: Surround by the self piece
                if (num_neighborhood_self + num_neighborhood_boun == MAXDIRECTION) {
                int check_flag = 0,
                check_eye_flag = num_neighborhood_boun;
                for (int d = 0 ; d < MAXDIRECTION; ++d) {
                    // Avoid fill self eye
                    if (NeighboorhoodState[d]==SELF && Liberties[d] > 1) {
                    check_eye_flag++;
                    }
                    // Check if there is one self component which has more than one liberty
                    if (NeighboorhoodState[d]==SELF && Liberties[d] > 1) {
                    check_flag = 1;
                    }
                }
                if (check_flag == 1 && check_eye_flag!=4) {
                    next_x = x;
                    next_y = y;
                }
                }
                // Case 2.2: Surround by opponent or both side's pieces.
                else if (num_neighborhood_oppo > 0) {
                int check_flag = 0;
                int eat_flag = 0;
                for (int d = 0 ; d < MAXDIRECTION; ++d) {
                    // Check if there is one self component which has more than one liberty
                    if (NeighboorhoodState[d]==SELF && Liberties[d] > 1) {
                    check_flag = 1;
                    }
                    // Check if there is one opponent's component which has exact one liberty
                    if (NeighboorhoodState[d]==OPPONENT && Liberties[d] == 1) {
                    eat_flag = 1;
                    }
                }
                if (check_flag == 1) {
                    next_x = x;
                    next_y = y;
                    if (eat_flag == 1) {
                    eat_move = 1;
                    }
                }
                else { // check_flag == 0
                    if (eat_flag == 1) {
                    next_x = x;
                    next_y = y;
                    eat_move = 1;
                    }
                }
                }
             }
         }
         bool repeat_move = 0;
		 if (next_x !=0 && next_y !=0) {
            if(game_length >= 10){
                // copy the current board to next board
                copy(&Board[0][0], &Board[0][0] + 11*11, &NextBoard[0][0]);
                // do the move
                // The move is a capture move and the board needs to be updated.
                if (eat_move == 1) {
                    update_board(NextBoard, next_x, next_y, turn);
                }
                else {
                    NextBoard[x][y] = turn;
                }
                // Check the history to avoid the repeat board (what the fuck?)
                bool repeat_flag = 1;
                if (NextBoard[next_x][next_y] != GameRecord[game_length-1][next_x][next_y]) {
                    repeat_flag = 0;
                }

                if (repeat_flag == 1 && eat_move == 1) {
                    repeat_move = 1;
                }
            }
		    if (repeat_move == 0) {
                // 3 digit zxy, z means eat or not, and put at (x, y)
                MoveList[legal_moves] = eat_move * 100 + next_x * 10 + y ;
                legal_moves++;
                if(random_bit == 1 && legal_moves == 1){
                    break;
                }
		    }
		 }
	    }
	    intersection_checked++;
    }
    return legal_moves;
}

int max_samples(int num_legal_moves, int first_N_round){
    int Num_max_samples = num_legal_moves*first_N_round;
    if(num_legal_moves>60)
        Num_max_samples += 10000;
    if(num_legal_moves<=60 && num_legal_moves>50)
        Num_max_samples += 15000;
    if(num_legal_moves<=50 && num_legal_moves>40)
        Num_max_samples += 18000;
    if(num_legal_moves<=40 && num_legal_moves>20)
        Num_max_samples += 20000;
    if(num_legal_moves<=20)
        Num_max_samples += 5000;
    return Num_max_samples;
}

float simulation(){
    return (rand()%10007) - 5000;
}

void Expansion(struct node *Root, int N , int samples){
    struct node *curr = NULL;
    struct node *prev = NULL;
    int i,j;
    // N: how many nodes to expand

    curr = (node *)malloc(sizeof(struct node));
    prev = (node *)malloc(sizeof(struct node));

    curr->firstchild = curr->nextsibling = NULL;
    curr->branch_id = 1;
    curr->games = samples;
    curr->score = 0;
    for(j=1 ; j<=samples ; j++){
        curr->score += simulation();
    }

    Root->firstchild = curr;
    prev = curr;

    for(i=2 ; i<=N ; i++){
        curr = (node *)malloc(sizeof(struct node));
        curr->firstchild = curr->nextsibling = NULL;
        curr->branch_id = i;
        curr->games = samples;
        curr->score = 0;
        for(j=1 ; j<=samples ; j++){
            curr->score += simulation();
        }

        prev->nextsibling = curr;
        prev = curr;
    }

    curr = Root->firstchild;
    while(curr != NULL){
        Root->games += samples;
        Root->score += curr->score;
        curr = curr->nextsibling;
    }
}

void simulate(int num_legal_moves, int first_N_round , int MoveList[HISTORYLENGTH],
              int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn,
              int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE],
              int scores_in_branch[NUMINTERSECTION], int games_in_branch[NUMINTERSECTION],
              double UCB[NUMINTERSECTION], int Num_max_samples){
    int i,j;
    int score;
    int sample = 0;
    int tmp;
    int move_id;
    int max_UCB_id = 0;
    double max_UCB = 0;

    int NextBoard[BOUNDARYSIZE][BOUNDARYSIZE]; // avoid to modify Board
    int New_game_length = 0;
    int UpdatedBoard[NUMINTERSECTION][BOUNDARYSIZE][BOUNDARYSIZE];
    int NewGameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE];
    int NewMoveList[HISTORYLENGTH]; // avoid to modify MoveList
    int turn_reset; // avoid to modify turn
    int next_legal_moves = 0; // avoid to modify num_legal_move

    if(game_length > 0 ){
        copy(&GameRecord[0][0][0], &GameRecord[0][0][0] + game_length*BOUNDARYSIZE*BOUNDARYSIZE,
             &NewGameRecord[0][0][0]);
    }

    for(sample=1 ; sample<=Num_max_samples+1 ; sample++){
        if(sample <= num_legal_moves*first_N_round)
            move_id = sample%num_legal_moves;
        // calculate max UCB
        if(sample > num_legal_moves*first_N_round){
            tmp = log(sample-1);
            for(i=0 ; i<num_legal_moves ; i++){
                UCB[i] = (double) scores_in_branch[i]/games_in_branch[i] + 1.414*sqrt(tmp/games_in_branch[i]);
                if(i==0){
                    max_UCB_id = 0;
                    max_UCB = UCB[0];
                }
                else if(UCB[i]>max_UCB){
                    max_UCB_id = i;
                    max_UCB = UCB[i];
                }

            }
            move_id = max_UCB_id;
            //fprintf(fp,"Best UCB Move id %d with UCB %.3f at length %d sample %d\n",move_id,UCB[move_id],game_length,sample);
        }
        if(sample == Num_max_samples+1)
            break;

        New_game_length = game_length;
        copy(&Board[0][0], &Board[0][0] + 11*11, &NextBoard[0][0]);
        turn_reset = turn;
        for(j=1 ; j<150 ; j++){
            // update
            if(j==1 && sample<=num_legal_moves){
                do_move(NextBoard, turn_reset, MoveList[move_id]);
                copy(&NextBoard[0][0], &NextBoard[0][0] + 11*11, &UpdatedBoard[move_id][0][0]);
                //fprintf(fp,"%d %d %d\n",(MoveList[move_id] % 100) / 10,MoveList[move_id] % 10, turn);
            }
            else if(j==1 && sample>num_legal_moves){
                copy(&UpdatedBoard[move_id][0][0], &UpdatedBoard[move_id][0][0] + 11*11, &NextBoard[0][0]);
                //fprintf(fp,"%d %d %d\n",(MoveList[move_id] % 100) / 10,MoveList[move_id] % 10, turn);
            }
            else{
                do_move(NextBoard, turn_reset, NewMoveList[0]);
            }
            // change turn
            if(turn_reset == BLACK){
                turn_reset = WHITE;
            }
            else{
                turn_reset = BLACK;
            }

            //copy(&NextBoard[0][0], &NextBoard[0][0] + BOUNDARYSIZE*BOUNDARYSIZE,
                 //&NewGameRecord[New_game_length][0][0]);
            //New_game_length++;
            next_legal_moves = gen_legal_move(NextBoard, turn_reset, New_game_length, NewGameRecord, NewMoveList, RANDOM);
            if(next_legal_moves == 0){
                score = final_score(NextBoard);
                // penalize lose points
                if(score <= DEFAULTKOMI){
                    scores_in_branch[move_id] = (float) scores_in_branch[move_id] + 25*(score-DEFAULTKOMI);
                    games_in_branch[move_id]++;
                }
                else {
                    scores_in_branch[move_id] = (float) scores_in_branch[move_id] + (score-DEFAULTKOMI);
                    games_in_branch[move_id]++;
                }
                break;
            }
            else if(next_legal_moves != 0)
                continue;
        }
    }
}
/*
 * This function randomly selects one move from the MoveList.
 * */
int rand_pick_move(int num_legal_moves, int MoveList[HISTORYLENGTH], int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    clock_t start_t, end_t;
    double time_taken;
    int i,j;
    int tmp;
    int sample = 0; // do samples
     // samples for each parent
    int first_N_round = 150; // UCB first n round
    int Num_max_samples = max_samples(num_legal_moves,first_N_round);
    // variables for statistic use
    int score;
    float max_score;
    double UCB[NUMINTERSECTION];
    int scores_in_branch[NUMINTERSECTION];
    int games_in_branch[NUMINTERSECTION];
    int max_move_id = 0;

    // NUMINTERSECTION

    // write something to file
    FILE *fp;
    fp = fopen("stdout.txt","a");
    setbuf(fp, NULL);
    srand(time(NULL));

    struct node *Root = NULL;

    Root = (node *)malloc(sizeof(struct node));

    Root->branch_id = 0; // root
    Root->games = 0;
    Root->score = 0;

    if (num_legal_moves == 0)
        return 0;
    else {
    // dingshi
    int move_id, move_eat, move_x, move_y;
    int star_x[9] = {3,3,3,5,5,5,7,7,7};
    int star_y[9] = {3,5,7,3,5,7,3,5,7};
    int star_valid[9] = {0,0,0,0,0,0,0,0,0};
    j = 0;
    for(i=0 ; i<num_legal_moves ; i++){
        move_id = MoveList[i];
        move_eat = move_id / 100;
        move_x = (move_id % 100) / 10;
        move_y = move_id % 10;

        if(move_x==3 | move_x==5 | move_x==7){
            if(move_y==3 | move_y==5 | move_y==7){
                star_valid[j] = move_eat*100 + move_x*10 + move_y;
                j++;
            }
        }
    }
    if(j>0 && game_length<7){
        tmp = rand()%j;
        return star_valid[tmp];
    }

    // simulate
	for(i=0 ; i<NUMINTERSECTION ; i++){
        scores_in_branch[i] = 0;
        games_in_branch[i] = 0;
        UCB[i] = 0;
	}

	start_t = clock();

	simulate(num_legal_moves, first_N_round, MoveList, Board, turn, game_length, GameRecord,
             scores_in_branch, games_in_branch, UCB, Num_max_samples);

	end_t = clock();
    time_taken = ((double)(end_t-start_t))/CLOCKS_PER_SEC;
    fprintf(fp,"Time: %.3f seconds at length %d, with legal moves: %d and samples: %d\n", time_taken, game_length, num_legal_moves, Num_max_samples);
	// statistics
	for(i=0 ; i<num_legal_moves ; i++){
        if(i == 0){
            max_score = UCB[i];
            max_move_id = i;
        }
        else if(UCB[i]>max_score){
            max_score = UCB[i];
            max_move_id = i;
        }
	}
	return MoveList[max_move_id];
    }
}
/*
 * This function update the Board with put 'turn' at (x,y)
 * where x = (move % 100) / 10 and y = move % 10.
 * Note this function will not check 'move' is legal or not.
 * */
void do_move(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int move) {
    int move_x = (move % 100) / 10;
    int move_y = move % 10;
    if (move<100) {
        Board[move_x][move_y] = turn;
    }
    else {
        update_board(Board, move_x, move_y, turn);
    }

}
/*
 * This function records the current game baord with current
 * game length "game_length"
 * */
void record(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE], int game_length) {
		for (int i = 0 ; i < BOUNDARYSIZE; ++i) {
		    for (int j = 0 ; j < BOUNDARYSIZE; ++j) {
			GameRecord[game_length][i][j] = Board[i][j];
		    }
		}
}
/*
 * This function randomly generate one legal move (x, y) with return value x*10+y,
 * if there is no legal move the function will return 0.
 * */
int genmove(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int turn, int time_limit, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    clock_t start_t, end_t, now_t;
    // record start time
    start_t = clock();
    // calculate the time bound
    end_t = start_t + CLOCKS_PER_SEC * time_limit;

    int MoveList[HISTORYLENGTH];
    int num_legal_moves = 0;
    int return_move = 0;

    num_legal_moves = gen_legal_move(Board, turn, game_length, GameRecord, MoveList, NONRANDOM);

    return_move = rand_pick_move(num_legal_moves, MoveList, Board, turn, game_length, GameRecord);

    do_move(Board, turn, return_move);

    return return_move % 100;
}
/*
 * This function counts the number of points remains in the board by Black's view
 * */
double final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    int black, white;
    black = white = 0;
    int is_black, is_white;
    for (int i = 1 ; i <= BOARDSIZE; ++i) {
	for (int j = 1; j <= BOARDSIZE; ++j) {
	    switch(Board[i][j]) {
		case EMPTY:
		    is_black = is_white = 0;
		    for(int d = 0 ; d < MAXDIRECTION; ++d) {
			if (Board[i+DirectionX[d]][j+DirectionY[d]] == BLACK) is_black = 1;
			if (Board[i+DirectionX[d]][j+DirectionY[d]] == WHITE) is_white = 1;
		    }
		    if (is_black + is_white == 1) {
			black += is_black;
			white += is_white;
		    }
		    break;
		case WHITE:
		    white++;
		    break;
		case BLACK:
		    black++;
		    break;
	    }
	}
    }
    return black - white;
}
/*
 * Following are commands for Go Text Protocol (GTP)
 *
 * */
const char *KnownCommands[]={
    "protocol_version",
    "name",
    "version",
    "known_command",
    "list_commands",
    "quit",
    "boardsize",
    "clear_board",
    "komi",
    "play",
    "genmove",
    "undo",
    "quit",
    "showboard",
    "final_score"
};

void gtp_final_score(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    double result;
    result = final_score(Board);
    result -= _komi;
    cout << "= ";
    if (result > 0.0) { // Black win
	cout << "B+" << result << endl << endl<< endl;;
    }
    if (result < 0.0) { // White win
	cout << "W+" << -result << endl << endl<< endl;;
    }
    else { // draw
	cout << "0" << endl << endl<< endl;;
    }
}
void gtp_undo(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    if (game_length!=0) {
	for (int i = 1; i <= BOARDSIZE; ++i) {
	    for (int j = 1; j <= BOARDSIZE; ++j) {
		Board[i][j] = GameRecord[game_length][i][j];
	    }
	}
    }
    cout << "= " << endl << endl;
}
void gtp_showboard(int Board[BOUNDARYSIZE][BOUNDARYSIZE]) {
    for (int i = 1; i <=BOARDSIZE; ++i) {
	cout << "#";
	cout <<10-i;
	for (int j = 1; j <=BOARDSIZE; ++j) {
	    switch(Board[i][j]) {
		case EMPTY: cout << " .";break;
		case BLACK: cout << " X";break;
		case WHITE: cout << " O";break;
	    }
	}
	cout << endl;
    }
    cout << "#  ";
    for (int i = 1; i <=BOARDSIZE; ++i)
	cout << LabelX[i] <<" ";
    cout << endl;
    cout << endl;

}
void gtp_protocol_version() {
    cout <<"= 2"<<endl<< endl;
}
void gtp_name() {
    cout <<"= TCG-randomGo99" << endl<< endl;
}
void gtp_version() {
    cout << "= 1.02" << endl << endl;
}
void gtp_list_commands(){
    cout <<"= ";
    for (int i = 0 ; i < NUMGTPCOMMANDS; ++i) {
	cout <<KnownCommands[i] << endl;
    }
    cout << endl;
}
void gtp_known_command(const char Input[]) {
    for (int i = 0 ; i < NUMGTPCOMMANDS; ++i) {
	if (strcmp(Input, KnownCommands[i])==0) {
	    cout << "= true" << endl<< endl;
	    return;
	}
    }
    cout << "= false" << endl<< endl;
}
void gtp_boardsize(int size) {
    if (size!=9) {
	cout << "? unacceptable size" << endl<< endl;
    }
    else {
	_board_size = size;
	cout << "= "<<endl<<endl;
    }
}
void gtp_clear_board(int Board[BOUNDARYSIZE][BOUNDARYSIZE], int NumCapture[]) {
    reset(Board);
    NumCapture[BLACK] = NumCapture[WHITE] = 0;
    cout << "= "<<endl<<endl;
}
void gtp_komi(double komi) {
    _komi = komi;
    cout << "= "<<endl<<endl;
}
void gtp_play(char Color[], char Move[], int Board[BOUNDARYSIZE][BOUNDARYSIZE], int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]) {
    int turn, move_i, move_j;
    if (Color[0] =='b' || Color[0] == 'B')
	turn = BLACK;
    else
	turn = WHITE;
    if (strcmp(Move, "PASS") == 0 || strcmp(Move, "pass")==0) {
	record(Board, GameRecord, game_length+1);
    }
    else {
	// [ABCDEFGHJ][1-9], there is no I in the index.
	Move[0] = toupper(Move[0]);
	move_j = Move[0]-'A'+1;
	if (move_j == 10) move_j = 9;
	move_i = 10-(Move[1]-'0');
	update_board(Board, move_i, move_j, turn);
	record(Board, GameRecord, game_length+1);
    }
    cout << "= "<<endl<<endl;
}
void gtp_genmove(int Board[BOUNDARYSIZE][BOUNDARYSIZE], char Color[], int time_limit, int game_length, int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]){
    int turn = (Color[0]=='b'||Color[0]=='B')?BLACK:WHITE;
    int move = genmove(Board, turn, time_limit, game_length, GameRecord);
    int move_i, move_j;
    record(Board, GameRecord, game_length+1);
    if (move==0) {
	cout << "= PASS" << endl<< endl<< endl;
    }
    else {
	move_i = (move%100)/10;
	move_j = (move%10);
//	cerr << "#turn("<<game_length<<"): (move, move_i,move_j)" << turn << ": " << move<< " " << move_i << " " << move_j << endl;
	cout << "= " << LabelX[move_j]<<10-move_i<<endl<< endl;
    }
}
/*
 * This main function is used of the gtp protocol
 * */
void gtp_main(int display) {
    char Input[COMMANDLENGTH]="";
    char Command[COMMANDLENGTH]="";
    char Parameter[COMMANDLENGTH]="";
    char Move[4]="";
    char Color[6]="";
    int ivalue;
    double dvalue;
    int Board[BOUNDARYSIZE][BOUNDARYSIZE]={{0}};
    int NumCapture[3]={0};// 1:Black, 2: White
    int time_limit = DEFAULTTIME;
    int GameRecord[MAXGAMELENGTH][BOUNDARYSIZE][BOUNDARYSIZE]={{{0}}};
    int game_length = 0;
    if (display==1) {
	gtp_list_commands();
	gtp_showboard(Board);
    }
    while (gets(Input) != 0) {
	sscanf(Input, "%s", Command);
	if (Command[0]== '#')
	    continue;

	if (strcmp(Command, "protocol_version")==0) {
	    gtp_protocol_version();
	}
	else if (strcmp(Command, "name")==0) {
	    gtp_name();
	}
	else if (strcmp(Command, "version")==0) {
	    gtp_version();
	}
	else if (strcmp(Command, "list_commands")==0) {
	    gtp_list_commands();
	}
	else if (strcmp(Command, "known_command")==0) {
	    sscanf(Input, "known_command %s", Parameter);
	    gtp_known_command(Parameter);
	}
	else if (strcmp(Command, "boardsize")==0) {
	    sscanf(Input, "boardsize %d", &ivalue);
	    gtp_boardsize(ivalue);
	}
	else if (strcmp(Command, "clear_board")==0) {
	    gtp_clear_board(Board, NumCapture);
	    game_length = 0;
	}
	else if (strcmp(Command, "komi")==0) {
	    sscanf(Input, "komi %lf", &dvalue);
	    gtp_komi(dvalue);
	}
	else if (strcmp(Command, "play")==0) {
	    sscanf(Input, "play %s %s", Color, Move);
	    gtp_play(Color, Move, Board, game_length, GameRecord);
	    game_length++;
	    if (display==1) {
		gtp_showboard(Board);
	    }
	}
	else if (strcmp(Command, "genmove")==0) {
	    sscanf(Input, "genmove %s", Color);
	    gtp_genmove(Board, Color, time_limit, game_length, GameRecord);
	    game_length++;
	    if (display==1) {
		gtp_showboard(Board);
	    }
	}
	else if (strcmp(Command, "quit")==0) {
	    break;
	}
	else if (strcmp(Command, "showboard")==0) {
	    gtp_showboard(Board);
	}
	else if (strcmp(Command, "undo")==0) {
	    game_length--;
	    gtp_undo(Board, game_length, GameRecord);
	    if (display==1) {
		gtp_showboard(Board);
	    }
	}
	else if (strcmp(Command, "final_score")==0) {
	    if (display==1) {
		gtp_showboard(Board);
	    }
	    gtp_final_score(Board);
	}
    }
}
int main(int argc, char* argv[]) {
    FILE *fp;
    fp = fopen("stdout.txt","w");
    fclose(fp);
//    int type = GTPVERSION;// 1: local version, 2: gtp version
    int type = GTPVERSION;// 1: local version, 2: gtp version
    int display = 0; // 1: display, 2 nodisplay
    if (argc > 1) {
	if (strcmp(argv[1], "-display")==0) {
	    display = 1;
	}
	if (strcmp(argv[1], "-nodisplay")==0) {
	    display = 0;
	}
    }
    gtp_main(display);
    return 0;
}
