#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

typedef struct matrice_caree{
    float **T;
} mat;
typedef struct vecteur_colonne{
    float **T;
} vect;

                                                          // REMARQUE :
            //j'ai préféré de travailler avec un vecteur colonne qui porte une seule colonne comme un tablaux de 2 Diemension
                                                //méme pour les structures matrices 


void temp_mat(float ** T,float r,int diag ,int ligne); //matrice temporaire pour stocker les valeurs de r 
float** alloc_mat(int n ,int m);//allocations dynamique des matrices
float** affect_mat(float ** T , int n ,int m);//l'initialisation des cases manuellement
float** rand_mat(float ** T , int n, int m );//l'initialisation des cases aleatoirement avec spécification de l'interval des valeurs genérer
float** dupliquer_mat(float**  A ,int n,int m);//formation d'une nouvelle matrice a partir d'une matrice donner comme argument pour appliquer des opérations arithmetique sans tocher la matrice principale
void aff_mat(float ** T ,int n, int m);//affichage d'une matrice quelconque
float** swap_row_mat(float** T , int diag , int ligne , int dim);//permutation des lignes d'une matrice lors de la rencontre d'une diagonal nul 
float** swap_row_vect(float** T , int diag , int ligne);//permutation des lignes d'un vecteur simultanement avec la fonction precedante .en but de resolution d'un system
float** rempl_ligne_mat(float** T ,int ligne ,int diag , int dim ,float g);//l'ajout d'une nouvelle valeur a une case spécifique d'une matrice 
float** rempl_ligne_vect(float** T ,int ligne ,int diag ,float g);//l'ajout d'une nouvelle valeur a une case spécifique d'un vecteur simultanement avec la fonction precedante .en but de resolution d'un system
void pivot_Gauss(float** T ,float** vect ,int n );//utilisation de la méthode de gauss pour avoir une matrice triangulaire supérieure 
float det(float ** T , int n );//calcul du determiant d'une matrice 
void init_mat(float** T,float** temp, int n);//l'initialisation de la matrice L triangulaire inférieure en utilisant la matrice temporaire ou on a stocker les val des cases inférieur au diagonal
void free_mat(float** T , int n );//la liberation de la mémoire qui a été alloer dynamiquement
void aff_prod_mat(float** A,float** U,float** L ,int n ,int p);//responsable de l'affichage de la decomposition LU
void aff_sys_X(float**A ,float**B,int n);//affiche une equation matriciale avec un vecteur colonne inconnu qui contient des X
void aff_sys_Y(float**L ,float**B,int n);//affiche une equation matriciale avec un vecteur colonne inconnu qui contient des Y (c'est justemet pour changer l'affichae si tu veux)
void resol_sys(float ** U ,float** X ,float **B ,int n );//resoudre un system qui contienne une matrice triangulaire supérieur et il ne retourne rien 
float** resol_ver(float ** U ,float** X ,float **B ,int n );//fait la méme chose que la fonction precedante mais cette fois avec un vecteur colonne comme type de retour
float** resol_inv(float ** U ,float** X ,float **B ,int n ); //resoudre un system qui contienne une matrice triangulaire inférieur et il retourne un vecteur colonne comme type de retour
void aff_sol (float**T ,int n );//affiche le vecteur colonne qui porte dans ces cases les valeurs solutions des inconnus Xs
float** init_base_canonique(float** T , int r , int n);//fait l'initialisation d'une base canonique qui porte un seul un en depandant de r donner comme argumant 
float** inverse_mat(float** U,float** L,int n );//spécialiser au calcule de l'inverse d'une matrice en utilisant la méthode numerique 
int verify_diag(float ** T,int n);//verifie les mineurs fendamentaux de la matrice pour qu'il admet une decomposition LU
float** dec_LU(float** T ,float** temp, int n );

float** alloc_mat(int n ,int m)
{
    int i ;
    float ** R =NULL ;
    R = (float**)malloc(sizeof(float*)*n);
    if(R==NULL)
    {
        printf("eurreure d'allocation 0\n");
        exit(33);
    }
    for (i=0;i<n;i++)
    {
        R[i] = (float*)malloc(sizeof(float)*m);
        if(R[i] == NULL)
        {
            printf("eurreure d'allocation %d \n",i+1);
            exit(34 + i);
        }
    }
    return R ;
}
float** affect_mat(float ** T , int n ,int m)
{
    int i , j ;
    for (i=0 ; i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            printf(" M[%d][%d] = ",i,j);
            scanf("%f",&T[i][j]);
        }
        printf("\n");
    }
    return T;
}
float** rand_mat(float ** T , int n,int m )
{
    int i , j , max,min ;
    if (m==1){
    printf("donner le plus grand element dans le vecteur :");
    scanf("%d",&max);
    printf("donner le plus petit element dans le vecteur :");
    scanf("%d",&min);
    }else{
    printf("donner le plus grand element dans la matrice :");
    scanf("%d",&max);
    printf("donner le plus petit element dans la matrice :");
    scanf("%d",&min);
    }
    for (i=0 ; i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            T[i][j]= min + rand()%(max - min + 1);
        }
    }
    return T;
}
void aff_mat(float ** T ,int n,int m)
{
    int i , j ;
    for (i=0;i<n;i++)
    {
        printf("( ");
        for(j=0;j<m;j++)
        {
            if(T[i][j]<0||T[i][j]>=10){printf(" %.2f %s",T[i][j],(j==m-1)?" ":"|");}
            else{printf(" %.3f %s",T[i][j],(j==m-1)?" ":"|");}
        }
        printf(")\n");
    }
}
float** swap_row_mat(float** T , int diag , int ligne , int dim)
{
    int i ;
    float temp ;
    for(i=0; i<dim;i++)
    {
        temp = T[diag][i];
        T[diag][i] = T[ligne][i];
        T[ligne][i] =temp;
    }

    return T ;
}
float** swap_row_vect(float** T , int diag , int ligne )
{
    float temp ;
    temp = T[diag][0];
    T[diag][0] = T[ligne][0];
    T[ligne][0] = temp ;
    return T ;
}
float** rempl_ligne_mat(float** T ,int ligne ,int diag , int dim ,float g)
{
    int i ; 
    for(i=0;i<dim;i++)
    {
        T[ligne][i] -=g*T[diag][i];
    }
    return T;
}
float** rempl_ligne_vect(float** T ,int ligne ,int diag ,float g)
{
    T[ligne][0] -= g * T[diag][0];
    return T ;
}
void pivot_Gauss(float** T ,float** vect , int n )
{
    int i , j ,k ; 
    float r ;
    for (i=0;i<n-1;i++)
    {
        if (T[i][i] == 0)
        {
            for (j=i+1;j<n;j++)
            {
                if (T[j][i]==0) continue;
                else{
                    T = swap_row_mat(T,i,j,n);
                    vect = swap_row_vect(vect,i,j);
                    printf("\n--------------------------------------------------permutation de lignes :-----------------------------------------------------------\n");
                    aff_sys_X(T,vect,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    break;
                }
            }
            if(j==n){
                printf("\nmatrice A est singuliere\n");
                exit(20);
            }
            else{
                i--;
                continue;
            }
        }
        for(k=i+1;k<n;k++)
        {
            r = T[k][i]/T[i][i];
            T = rempl_ligne_mat(T,k,i,n,r);
            vect = rempl_ligne_vect (vect,k,i,r);
            if(k<n-1)
            {           
                printf("\n------------------------------------------------- processus en cours :--------------------------------------------------------------\n\n");
                aff_sys_X(T,vect,n);
                printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
            }
        }
    }
    printf("\n------------------------------------------------- processus terminer :--------------------------------------------------------------\n\n");
    aff_sys_X(T,vect,n);
    printf("-------------------------------------------------------------------------------------------------------------------------------------\n\n");
}

float det(float ** T , int n )
{
    int i ; 
    float det=1;
    for (i=0;i<n;i++)
    {
        det *= T[i][i];
    }
    return det ;
}
void init_mat(float** T,float** temp, int n)
{
    int i , j ;
    for (i=0;i<n;i++){
        for(j=0;j<n;j++){
            if((i<=j)){T[i][j]=((i<j)?0:1);}
            else {T[i][j]=temp[i][j];}
        }
    }
}
void temp_mat(float ** T,float r,int diag ,int ligne)
{
    T[ligne][diag] = r;
}
void free_mat(float** T , int n )
{
    int i ;
    for (i=0;i<n;i++)
    {
        free(T[i]);
    }
    free(T);
}
void aff_prod_mat(float** A,float** U,float** L ,int n ,int p)
{
    int i ,j ,k,m; 
    for (i=0;i<n;i++)
    {
        printf("(");
        for(j=0;j<p;j++)
        {
            if(A[i][j]<0||A[i][j]>=10){printf(" %.2f %s",A[i][j],(j==p-1)?"":"|");}
            else{printf(" %.3f %s",A[i][j],(j==p-1)?"":"|");}
        }
        printf(") %s ",(i==((int)n/2))?"=":" ");
        printf("(");
        for(j=0;j<p;j++)
        {
           if(L[i][j]<0||L[i][j]>=10){printf(" %.2f %s",L[i][j],(j==p-1)?"":"|");}
           else{printf(" %.3f %s",L[i][j],(j==p-1)?"":"|");}
        }
        printf(") %s ",(i==((int)n/2))?"*":" ");
        printf("(");
        for(j=0;j<p;j++)
        {
            if(U[i][j]<0||U[i][j]>=10){printf(" %.2f %s",U[i][j],(j==p-1)?"":"|");}
            else{printf(" %.3f %s",U[i][j],(j==p-1)?"":"|");}
        }
        printf(")\n");

    }
}
float** dupliquer_mat(float** A ,int n,int m)
{
    int i ,j ;
    float** U = alloc_mat(n,m);
    for (i=0;i<n;i++)
    {
        for(j=0;j<m;j++)
        {
            U[i][j]=A[i][j];
        }
    }
    return U ; 
}
void aff_sys_X(float**A ,float**B,int n)
{
    int i , j ;
    for (i=0;i<n;i++)
    {
        printf("(");
        for(j=0;j<n;j++)
        {
            if(A[i][j]<0||A[i][j]>=10){printf(" %.2f %s",A[i][j],(j==n-1)?"":"|");}
            else{printf(" %.3f %s",A[i][j],(j==n-1)?"":"|");}
        }
        printf(") %s ",(i==((int)n/2))?"*":" ");
        printf("( x%d",i+1);
        printf(" ) %s ",(i==((int)n/2))?"=":" ");
        if(B[i][0]<0||B[i][0]>10)
        {
            printf("( %.2f",B[i][0]);
        }
        else{
            printf("( %.3f",B[i][0]);
            }
        printf(" )\n");
    }
}
void aff_sys_Y(float**L ,float**B,int n)
{
    int i , j ;
    for (i=0;i<n;i++)
    {
        printf("(");
        for(j=0;j<n;j++)
        {
            if(L[i][j]<0||L[i][j]>=10){printf(" %.2f %s",L[i][j],(j==n-1)?" ":"|");}
            else{printf(" %.3f %s",L[i][j],(j==n-1)?" ":"|");}
        }
        printf(" ) %s ",(i==((int)n/2))?"*":" ");
        printf("( y%d",i+1);
        printf(" ) %s ",(i==((int)n/2))?"=":" ");
        if(B[i][0]<0||B[i][0]>10)
        {
            printf("( %.2f",B[i][0]);
        }
        else{
            printf("( %.3f",B[i][0]);
            }
        printf(" )\n");
    }
}
void resol_sys(float ** U ,float** X ,float **B ,int n )
{
    int i , j ;
    for (i=n-1;i>=0;i--)
    {
        X[i][0] = B[i][0];
        for(j=i+1;j<n;j++)
        {
            X[i][0] -= U[i][j] * X[j][0];
        }
        X[i][0] /= U[i][i];
    }
}
float** resol_ver(float ** U ,float** X ,float **B ,int n )
{
    int i , j ;
    for (i=n-1;i>=0;i--)
    {
        X[i][0] = B[i][0];
        for(j=i+1;j<n;j++)
        {
            X[i][0] -= U[i][j] * X[j][0];
        }
        X[i][0] /= U[i][i];
    }
    return X ;
}
float** resol_inv(float ** U ,float** X ,float **B ,int n )
{
    int i , j ;
    for (i=0;i<n;i++)
    {
        X[i][0] = B[i][0];
        for(j=0;j<i;j++)
        {
            X[i][0] -= U[i][j] * X[j][0];
        }
        X[i][0] /= U[i][i];
    }
    return X ;
}
void aff_sol (float**T ,int n  )
{
    int i ;
    for(i=0;i<n;i++)
    {
        printf("\t\t( x%d",i+1);
        printf(" ) = ");
        printf("(");
        if(T[i][0]<0||T[i][0]>=10){printf(" %.2f",T[i][0]);}
            else{printf(" %.3f",T[i][0]);}
        printf(" )\n");
    }
}
float** init_base_canonique(float** T , int r , int n)
{
    int i ;
    for (i=0;i<n;i++)
    {
        T[i][0]=(i==r)?1:0;
    }
    return T ; 
}
float** inverse_mat(float** U,float** L,int n )
{
    float** inverse = alloc_mat(n, n);
    int i ,j;
    for (i=0;i<n;i++)
    {
        float** e = alloc_mat(n, 1);
        float** Y = alloc_mat(n, 1); 
        float** X = alloc_mat(n, 1); 
        float** colo = alloc_mat(n, 1); 
        for (j=0;j<n;j++)
        {
            e = init_base_canonique(e,i,n);
            X = resol_inv(L,X,e,n);
            Y = dupliquer_mat(X,n,1);
            colo = resol_ver(U,colo,Y,n);
            inverse[j][i] = colo[j][0];
        }
        free_mat(e, n); // Libération de la mémoire pour e
        free_mat(Y, n); // Libération de la mémoire pour Y
        free_mat(X, n); // Libération de la mémoire pour X
        free_mat(colo, n); // Libération de la mémoire pour colo
    }
    return inverse ;
}

int verify_diag(float ** T,int n)
{
    int i ;
    for(i=0;i<n;i++)
    {
        if (T[i][i] == 0) return 0 ;
    }
    return 1 ;
}
float** dec_LU(float** T ,float** temp, int n )
{
    int i , j ,k ; 
    float r ;
    for (i=0;i<n-1;i++)
    {
        if (T[i][i] == 0)
        {
            for (j=i+1;j<n;j++)
            {
                if (T[j][i]==0) continue;
                else{
                    T = swap_row_mat(T,i,j,n);
                    break;
                }
            }
            if(j==n){
                printf("\nmatrice A est singuliere\n");
                exit(20);
            }
            else{
                i--;
                continue;
            }
        }
        for(k=i+1;k<n;k++)
        {
            r = T[k][i]/T[i][i];
            temp_mat(temp,r,i,k);
            T = rempl_ligne_mat(T,k,i,n,r);
        }
    }
    return T ;
}

int main ()
{
    //initialisation d'un generateur des nombres pseudo aléatoires au temps NULL 
    srand(time(NULL));
    int n,p=0,a=0,m=1,f=0;
    char e ,w,c;
    mat A,L,U,F,inverse,temp ;
    vect B,S,sol ;

    printf("\n          Bonjour , Ce programme a ete realiser par EL FARISSI Oussama Etudiant Premiere annee GI a EST Agadir .il permet de calculer:\n");
    printf ("___________________Determinant A /Pivot de Gauss/ Resolution d'un System /decomposition LU /calculer Inverse d'une matrice_________________\n");
    do{
        if (f>0)//il a au moins modifier une seule fois
        {
            printf("------------------------------------------------------------------------------------------------------------------------------------\n");
            printf ("\n\tla modification numero %d \n\n", f);
            printf("------------------------------------------------------------------------------------------------------------------------------------\n");
        }
        //l'intialisation de la modification d'eviter de tomber dans une boucle infinie
        m=0;
        do
            {
            printf ("\nentrer les dimension de la matrice carree A :");
            if(scanf("%d",&n) != 1){// Vérification de la saisie ,car scanf return != 1 en cas d'une saisie non valide(non entier)
                printf("Erreur : Veuillez entrer un entier valide.\n");
                scanf("%c", &c);// Lire et supprimer le caractère non entier du flux d'entrée
                n=0;// Réinitialiser n pour rester dans la boucle
            }
        }while(n<=1);//les dimensions doit étre strictement supérieur a 1 c-à-d au moins n == 2
        A.T = alloc_mat(n,n);
        B.T = alloc_mat(n,1);
        do
        {
            p=0;//l'initialisation d'erreur pour éviter de tomber dans une boucle infinie
            printf("voullez vous remplir la matrice A Manuelement ou Alleatoirement (M/A):");
            scanf(" %c",&e);
            switch(e)
            {
                case 'M' : 
                    affect_mat(A.T,n,n);
                break;
                case 'A' :
                    rand_mat(A.T,n,n);
                break;
                default:
                    printf("\adesole il y a une erreur de saisie vous voullez ressayer....\n");
                    p=1;//spécfieur d'erreur
            }
        }while(p==1);//la repetition de la boucle do while en cas d'une erreur de saisie
        do
        {
            printf("voullez vous remplir le vecteur B Manuelement ou Alleatoirement (M/A):");
            scanf(" %c",&e);
            switch(e)
            {
                case 'M' : 
                    affect_mat(B.T,n,1);
                break;
                case 'A' :
                    rand_mat(B.T,n,1);
                break;
                default:
                    printf("\adesole il y a une erreur de saisie vous voullez ressayer....\n");
                    p=1;//spécfieur d'erreur
            }
        }while(p==1);//la repetition de la boucle do while en cas d'une erreur de saisie
        // printf("\n-------------------------------------------------------------- vecteur B :----------------------------------------------------------\n\n");
        // aff_mat(B.T,n,1);
        // printf("------------------------------------------------------------------------------------------------------------------------------------\n");
        temp.T =alloc_mat(n,n);
        L.T = alloc_mat(n,n);
        sol.T = alloc_mat(n,1);
        S.T =dupliquer_mat(B.T,n,1);
        F.T =dupliquer_mat(A.T,n,n);
        U.T =dupliquer_mat(A.T,n,n);
        U.T =dec_LU(U.T,temp.T,n);
        init_mat(L.T,temp.T,n);
        inverse.T = inverse_mat(U.T,L.T,n);
        do{
            m=0;//l'intialisation de la modification pour éviter de tomber dans une boucle infinie
            printf("\n-------------------------------------voullez selectionner une tache a executer :----------------------------------------------------\n");
            printf("(0=Quitter) (1=Determinant) (2=Pivot de Gauss) (3=Resolution system) (4=decomposition LU) (5=Inverse d'une matrice) (6=modification)\n");
            printf("------------------------------------------------------------------------------------------------------------------------------------\n");
            printf("entrer votre choix ici :");
            scanf("%s",&w);
            switch(w)
            {
                case '0' :
                    a=1;// l'utilisateur veut quitter
                break;
                case '1' :
                    printf("\n-------------------------------------------------------------- matrice A :----------------------------------------------------------\n\n");
                    aff_mat(A.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    printf("---------------------------------------------------le determinant de la matrice A est : --------------------------------------------\n");
                    printf("Det(A) = %.2f\n",det(U.T,n));
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                break;
                case '2' :
                    printf("\n------------------------------------------------------------ Pivot de Gauss :-------------------------------------------------------\n\n");
                    pivot_Gauss(F.T,S.T,n);
                    // aff_mat(U.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                break ;
                case '3' :
                    printf("\n---------------------------------------------------- Resolution de System AX = B :---------------------------------------------------\n\n");
                    aff_sys_X(A.T,B.T,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    resol_sys(U.T,sol.T,S.T,n);
                    aff_sol(sol.T,n);
                break;
                case '4' :
                    if (det(U.T,n)==0 ) {
                        printf ("\n-----------------------------------matrice A est non inversible et son determinant egal a zero  ----------------------------------\n\n");
                        printf ("-------------------vous voullez selectionner une autre tache a executer ou bien d'entrer une autre matrice inversible -------------\n\n");
                    }
                    else if(verify_diag(A.T,n) == 0 ){
                        printf ("\n------------------------------les mineurs fondamentaux ne sont pas tous differents de zero --------------------------------------\n\n");
                        printf ("-------------------vous voullez selectionner une autre tache a executer ou bien d'entrer une autre matrice inversible -------------\n\n");
                    }
                    else{
                    printf("\n--------------------------------------------la matrice A admet un decomposition LU :-------------------------------------------------\n\n");
                    aff_prod_mat(A.T,U.T,L.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    printf("\n-------------------------------------------------------------- matrice L :----------------------------------------------------------\n\n");
                    aff_mat(L.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    printf("\n-------------------------------------------------------------- matrice U :----------------------------------------------------------\n\n");
                    aff_mat(U.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    }

                break;
                case '5' :
                    printf("\n------------------------------------------------------------ matrice A :------------------------------------------------------------\n\n");
                    aff_mat(A.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    printf("\n-------------------------------------------------------l'inverse de A est :-------------------------------------------------------\n\n");
                    aff_mat(inverse.T,n,n);
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                break;
                case '6' :
                    m=1;//l'utilisateur veut modifier 
                    f++;//compteur de nombre de modifications  
                break;
                default:
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
                    printf("\n\tdesole il y a une eureur de saisie voullez ressayer .......\n");
                    printf("------------------------------------------------------------------------------------------------------------------------------------\n");
            }
        //a==0 => il ne veut pas quitter et m==0 => il ne veut pas modifier 
        }while(a==0 && m==0);
    //a==0 => il ne veut pas quitter ou m==1 => il veut modifier 
    }while(a==0||m==1);
    printf("\n----------------------------------------------------merci pour votre attention ---------------------------------------------------\n");
    
    //liberation de la mémoire
    free_mat(A.T,n);
    free_mat(L.T,n);
    free_mat(U.T,n);
    free_mat(temp.T,n);
    free_mat(B.T,1);
    free_mat(inverse.T,n);
    free_mat(sol.T,1);
    free_mat(F.T,n);
    return 0;
}