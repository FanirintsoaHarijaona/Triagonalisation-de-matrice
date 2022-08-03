#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
//bibliothèque inclue pour tester le type récupéré dans le fichier
#include <typeinfo>
using namespace std;

vector<vector<float>> triangulariser(vector<vector<float>>& matrice,vector<float> &b,int& dim);
void afficherMatrice(vector<vector<float>>& matrice,int& dim);
vector<float> pivot(vector<vector<float>>& matrice, int colonne,int dim);
void echangeLignes(vector<float>&matricePivotInitial,vector<float>&matricePivotMax,int &dim);
void initialiserMatrice(vector<vector<float>> &matrice,vector<float> &second,int &dim);
void resolution(vector<vector<float>> &matrice,vector<float> &second,int &dim);

int main(){
    cout<<"Triangularisation par la méthode de Gauss de la matrice\n";
    vector<vector<float>> matrice;
    vector<float> second;
    //initialisation des données:matrice,dimension
    int dim = 0;
    initialiserMatrice(matrice,second,dim);
    cout<<"La matrice:"<<endl;
    afficherMatrice(matrice,dim);
    cout<<endl;
    vector<vector<float>> result = triangulariser(matrice,second,dim);
    cout<<"Matrice triangularisée :"<<endl;
    afficherMatrice(result, dim);
    cout<<"Résolution de l'équation"<<endl;
    resolution(matrice,second,dim);
    return 0;
}

void initialiserMatrice(vector<vector<float>> &matrice,vector<float> &second,int &dim){
//initialisation d'une variable afin d'ouvrir le fichier
    ifstream fichier("data.txt");
    if(fichier){
//la première ligne du fichier est déstinéé à être la dimension de la matrice
        fichier>>dim;
//Tant que la variable i sera inférieure à la dimension insérée ci-dessus
//nous considérerons les lignes en dessous de la première comme étant les
//coéfficients de la matrice
        for(int i=0;i<dim;i++){
            vector<float> ligne;
            for(int j=0;j<dim;j++){
                float temp(0);
                fichier>>temp;
                ligne.push_back(temp);
            }
            matrice.push_back(ligne);
        }
//le reste du fichier sera la matrice de B du système d'équation AX = B où A est 
//la matrice à triangulariser
        for(int i=0;i<dim;i++){
            float temp = 0;
            fichier>>temp;
            second.push_back(temp);
        }
    }
//En cas d'erreur d'ouverture du fichier, on affiche un message d'erreur
else{
        cout<<"Erreur lors de chargement du fichier"<<endl;
    }
}

vector<vector<float>> triangulariser(vector<vector<float>>& matrice,vector<float>&b,int& dim){ 
    if(matrice[0].size()!=dim){
        cout<<"Matrice non carrée donc impossible à triangulariser"<<endl;
    }
     for(int k = 0;k< dim-1; k++){
//Si la ligne du pivot ne correspond à la diagonale ,
//on intervertit la ligne du pivot avec la ligne de l'élément de la diagonale
         int lignePivot =(int)pivot(matrice,k,dim)[1];
         cout<<"lignePivot:"<<lignePivot<<endl;
         if(lignePivot!=k){
            echangeLignes(matrice[k],matrice[lignePivot-1],dim);
//variable temporaire pour stocker la valeur de la matrice second membre
//Quand on intervertit les lignes de la matrice A, on intevertit aussi celui
//du seconde membre de l'équation
            float tmp = b[k];
            cout<<"tmp:"<<tmp;
            cout<<"Pivot:"<<b[lignePivot-1];
            b[k] = b[lignePivot-1];
            b[lignePivot-1] = tmp;
         }
/***après le choix du pivot et interversion,on calcule les autres coefficients
    les coefficients en dessous du pivot seront nuls
    tandis que ceux des autres sont calculées à partir du formule
    matrice[i][j] = matrice[i][j] - (matrice[i][k]/matrice[k][k])*matrice[k][j];
    (matrice[i][k]/matrice[k][k]) est le coefficient qui permet d'annuler les 
    coefficients en dessous du pivot
    ***/
        for (int i= k+1; i<dim ;i++ ){
            for (int j= k+1 ;j<dim;j++){
                matrice[i][j] = matrice[i][j] - (matrice[i][k]/matrice[k][k])*matrice[k][j];
            }
            b[i] = b[i] - (matrice[i][k]/matrice[k][k])*b[k];
            matrice[i][k] = 0;
            }
        }
//affichage du second membre une fois que la triangularisation est terminée
        cout<<"triangularisation achevée.Le second membre de l'équation est comme suit("<<endl;
        for(int i=0;i<b.size();i++){
        cout<<b[i]<<"\t";
    }
    cout<<")"<<endl<<endl<<endl;
    return matrice;
     }

void resolution(vector<vector<float>> &matrice,vector<float> &second,int &dim){
//fonction résolution du l'équation AX = B
    
//vecteur qui contient les inconnues du système
    vector<float> x ;
//initialisation des inconnues
    for(int l=0; l<dim;l++){
        x.push_back(0);
    }
    float somme;

    for(int i = dim-1; i>=0; i-- ){
        somme = 0;
        for(int j =i+1; j<dim; j++){
//produit de la matrice A.X
            somme+=(x[j]*matrice[i][j]);
        }
//On obtient le résultat grâce à b[i]-somme/elt de la diagonale de A
        x[i] = (second[i] - somme )/matrice[i][i];
    }
    cout<<"La résolution de l'équation AX = B(A =matrice , B=second,x:inconnu) :"<<endl;
    cout<<"x:(";
//affichage des solutions du système d'équations AX =B
    for(int i=0;i<x.size();i++){
        cout<<x[i]<<"\t";
    }
    cout<<")"<<endl;
}

vector<float> pivot(vector<vector<float>>& matrice, int ligne,int dim){
//fonction qui choisit le pivot
    vector <float> result;
    float max = matrice[0][ligne];
    float lignePivot = 0;
//Les lignes possédant déja un pivot ne sont plus analysées donc reste inchangées
//après l'interversion
    for (float i=ligne;i<dim;i++){
        if(max <fabs(matrice[i][ligne])){
            max = matrice[ligne][i];
//ligne du plu grand dans la colonnes
            lignePivot = i+1;
            }
    }
    result.push_back(max);
    result.push_back(lignePivot);
    return result;
}

void echangeLignes(vector<float>&matricePivotInitial,vector<float>&matricePivotMax,int &dim){
//echange des éléments de 2 lignes
    for(int i = 0; i<dim;i++){
        float tmp = 0;
        tmp =matricePivotInitial[i];
        matricePivotInitial[i]= matricePivotMax[i];
        matricePivotMax[i] = tmp;
    }   
}

//fonction qui affiche une matrice à 2 dimensions
void afficherMatrice(vector<vector<float>>& matrice,int& dim){
    for (int i=0;i<dim;i++){
       for (int j=0; j<dim;j++){
            cout<<matrice[i][j]<<"\t";
        }
       cout<<endl;
    }
}