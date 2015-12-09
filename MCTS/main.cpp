#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

using namespace std;

struct node {
  struct node *firstchild;
  struct node *nextsibling;
  int branch_id;
  int games;
  float score;
};

float simulation(){
    return (rand()%10007) - 5000;
}

struct node *find_path(struct node *Root){
    struct node *curr = NULL;
    struct node *best = NULL;
    struct node *path = NULL;
    float score = 0;
    float max_average_score = 0;
    int max_child_id = 0;

    //path = (node *)malloc(sizeof(struct node));
    best = (node *)malloc(sizeof(struct node));

    if(Root->firstchild==NULL){
        return NULL;
    }

    curr = Root->firstchild;
    score = (float) curr->score/curr->games;
    max_average_score = score;
    max_child_id = curr->branch_id;
    best = curr;

    while(curr != NULL){
        printf("...%d %d %.f\n",curr->branch_id, curr->games, curr->score);
        score = (float) curr->score/curr->games;
        if(score > max_average_score){
            max_average_score = score;
            max_child_id = curr->branch_id;
            best = curr;
        }
        curr = curr->nextsibling;
    }
    printf("Best: %d %d %.f\n",best->branch_id, best->games, best->score);
    best->firstchild = find_path(best);

    return best;
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

int main()
{
    int i,j;
    int N_total = 3;
    int Node_expansion = 5;
    int samples = 5;

    float update_score = 0;
    int update_games = 0;

    struct node *Root = NULL;
    struct node *curr = NULL;
    struct node *path = NULL;
    struct node *expansion_node = NULL;

    Root = (node *)malloc(sizeof(struct node));
    curr = (node *)malloc(sizeof(struct node));

    // Obtain an initial tree
    Root->branch_id = 0; // root
    Root->games = 0;
    Root->score = 0;

    Expansion(Root,Node_expansion,samples);

    printf("Root: %d Games: %d Score: %.f\n",Root->branch_id, Root->games, Root->score);

    for(i=1 ; i<=N_total ; i++){
        // selection
        printf("Selection:\n");
        path = find_path(Root);
        curr = path;
        while(curr != NULL){
            expansion_node = curr;
            printf("%d %d %.f\n",expansion_node->branch_id, expansion_node->games, expansion_node->score);
            curr = curr->firstchild;
        }
        // Expansion and Simulation
        printf("Expansion and Simulation:\n");
        update_score = expansion_node->score;
        update_games = expansion_node->games;
        Node_expansion = 5;
        samples = 2;
        Expansion(expansion_node,Node_expansion,samples);
        update_score = expansion_node->score - update_score;
        update_games = expansion_node->games - update_games;
        printf("%.f %d\n",update_score,update_games);
        // back propagation
        printf("Back Propagation:\n");
        curr = path;
        while(curr != expansion_node){
            curr->games += update_games;
            curr->score += update_score;
            printf("Updated: %d %d %.f\n",curr->branch_id, curr->games, curr->score);
            curr = curr->firstchild;
        }
    }
    return 0;
}



