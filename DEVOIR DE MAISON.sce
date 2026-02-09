//PROGRAMMATION DES METHODES CLASSIQUES
//DECLARATION DES FONCTIONS
//FONCTION POUR INVERSER UNE MATRICE
function [inversedematrice] = inverse_matrice(matrice)
        n=size(matrice,1)
    det_matrice = det(matrice);
    
    if det_matrice == 0
        disp("La matrice est singulière, son inverse n existe pas.")
    else
        inversedematrice = zeros(n,n);
        
        for i = 1:n
            for j = 1:n
                A = matrice;
                A(j,:) = [];
                A(:,i) = [];
                inversedematrice(i,j) = (-1)^(i+j) * det(A) / det_matrice;
            end
        end
    end
endfunction
funcprot(0)
function [x,r_kj,erreurJ,rspec]=Jacobi(A,b,c,epsi,n_max_iter)
        // decomposition de la matrice A = D-E-F
        n = size(A,1);
        // Initialisation de la matrice D avec des zéros
D = zeros(size(A,1), size(A,2))

// Remplissage de la diagonale de D avec les éléments diagonaux de A
for i=1:size(A,1)
    D(i,i) = A(i,i)
end

// Initialisation de la matrice E avec des zéros
E = zeros(size(A,1), size(A,2))

// Remplissage de la partie triangulaire inférieure de E
for i=1:size(A,1)
    for j=i:size(A,2)
        E(i,j) = -A(i,j)
    end
    E(i,i) = 0
end

// Initialisation de la matrice F avec des zéros
F = zeros(size(A,1), size(A,2))

// Remplissage de la partie triangulaire inférieure et diagonale de F
for i=1:size(A,1)
    for j=1:size(A,2)
        if j <= i
            F(i,j) = -A(i,j)
        end
        F(i,i) = 0
    end
end

 //CALCULONS L'INVERSE DE D SANS UTILISER LA FONCTION INV
for i=1:n
    inverseD(i,i)=1/D(i,i)
end
//DEFINITION DES MATRICES M,N et J
inverseM=inverseD
N=E+F
J=inverseM*N
// CALCULONS LE RAYON SPECTRALE J
rspec=max(abs(spec(J)))
if rspec<1 then
   printf("Il y a convergence \n")
end
if rspec >=1 then
   printf("Il n y a pas de convergence \n ")
            break;
end
// initialisation
n_iter = 1
//algorithme de Jacobi
x_kj=c

for k=1:n_max_iter
    x_kj=[x_kj c]
    c=inverseM*(N*c + b)
    if max(abs(A*c-b))<epsi
        n_iter = k
        break
     end
end
mat=x_kj(:,2:$)
        //VECTEUR RESIDU
for i=1:n_iter
    r_k(1,i)=(norm((A*(mat(:,i))-b),2))./(norm(b,2))
end
        // ERREUR RELATIVE
x=A\b
for i=1:n_iter
    erreurJ(1,i)=(norm(((mat(:,i))-x),2))./(norm(b,2))
end
        
        // Affichage des resultats 
        printf("Le rayon spectrale recherché  est : %f",rspec);
        disp("La solution est la suivante :",x);
        disp("La suite des résidus relatifs est : ",r_kj);
        disp("la suite des erreurs relatives est : ",erreurJ);
        disp("Nombre d iterration",n_iter);
endfunction

//FONCTION POUR INVERSER UNE MATRICE
function inversedematrice = inverse_matrice(matrice)
        n=size(matrice,1)
    det_matrice = det(matrice);
    
    if det_matrice == 0
        disp("La matrice est singulière, son inverse n existe pas.")
    else
        inversedematrice = zeros(n,n);
        
        for i = 1:n
            for j = 1:n
                A = matrice;
                A(j,:) = [];
                A(:,i) = [];
                inversedematrice(i,j) = (-1)^(i+j) * det(A) / det_matrice;
            end
        end
    end
endfunction

//DECLARATION DE LA FONCTION DE GAUSS-SEIDEL

    // fonction Gauss-Seidel
function [x,r_k,erreurG,rspec]=GaussSeidel(A,b,c,epsi,n_max_iter)
    // decomposition de A = D-E-F
    n = size(A,1);
    D= diag(diag(A))
    E=-triu(A)+D
    F=-tril(A)+D
    // definition des matrices M,N et J
    M= D-E
    N= F
    J=inverse_matrice(M)*N
    rspec=max(abs(spec(J)))
    if rspec<1 then
        printf("Il y a convergence \n")
    end
    
    if rspec >=1 then
        printf("Il n y a pas de convergence \n")
        abort;;
    end
    // initialisation
    n_iter = 0 
    //algorithme de Gauss
    // intialisation
    n_iter=0
    info=0
    x_k=c
    for k=1:n_max_iter
        x_k=[x_k c]
        c=inverse_matrice(M)*(N*c + b)
if max(abs(A*c-b))<epsi
   n_iter = k
   info=1
   break
end
end
if info==0 then
   printf("Vecteur n a pas convergé pour %d iterrations ",n_max_iter)
   abort;;
end
    mat=x_k(:,2:$)
    // vecteur residu
for i=1:n_iter
    r_k(1,i)=(norm((A*(mat(:,i))-b),2))./(norm(b,2))
end
    // vecteur e_kur
    x=A\b
for i=1:n_iter
     erreurG(1,i)=(norm(((mat(:,i))-x),2))./(norm(b,2));
end
    // AFFICHAGE DES RESULTATS
    printf("Le rayon spectrale est : %f",rspec);
    disp("La solution est :",x);
    disp("La suite des résidus relatifs est : ",r_k);
    disp("la suite des erreurs relatives est : ",erreurG);
    disp("Nombre d iterration",n_iter);
endfunction




// DECLARATION DE LA FONCTION DE L'ALGORITHME DE LA METHODE DE RELAXATION
 function [x,r_kR,erreurR,rspec]=Relaxation(A,b,c,epsi,omega,n_max_iter)
        // decomposition de A = D-E-F
        n = size(A,1);
        D= diag(diag(A))
        E=-triu(A)+D
        F=-tril(A)+D
        // deifinition des matrices M,N et J
        M=(1/omega)*D -E
        N=((1-omega)/omega)*D+F
        R_omega=inverse_matrice(omega*M)*(omega*N)
        rspec=max(abs(spec(R_omega)))
        if rspec<1 then
            printf("Il y a convergence")
        end
        if rspec >=1 then
            printf("Il n y a pas de convergence")
            abort;
        end
        //algorithme de Relaxation
        x_kR=c
        for k=1:n_max_iter
            x_kR=[x_kR c]
            c=inv(omega*M)*((omega*N)*c + omega*b)
            if max(abs(A*c-b))<epsi
                n_iter = k
                break
            end
        end
        mat=x_kR(:,2:$)
        // vecteur residu
    for i=1:n_iter
        r_kR(1,i)=(norm((A*(mat(:,i))-b),2))./(norm(b,2))
    end
    // vecteur e_kur
    x=A\b
    for i=1:n_iter
        erreurR(1,i)=(norm(((mat(:,i))-x),2))./(norm(b,2))
    end
    // Affichage
    printf("Le rayon spectrale est : %f",rspec)
    disp("La solution est :",x)
    disp("La suite des résidus relatifs est : ",r_kR)
    disp("la suite des erreurs relatives est : ",erreurR)
    disp("Nombre d iterration",n_iter)
endfunction



//fonction tridiaggAb
funcprot(0)
function [A,b,c]= tridiagAb(n,alpha,diago,g_amma)
A = zeros(n, n);

for i=1:n
    for j=1:n
        if i == j
            A(i,j) = diago;
        elseif i == j+1
            A(i,j) = g_amma;
        elseif i+1 == j
            A(i,j) = alpha;
        end
    end
end
b=zeros(n,1)
for i=1:n
    b(i) = sin(4*%pi*i/(n+1));
end
c=zeros(n,1)
endfunction
//PROGRAMME PRINCIAPAL ET APPLICATIONS DES FONCTIONS
    printf("[======Bienvenue dans le programme de Anon et Juane========]\n")
    printf("[-----------------CHOISISSEZ-----------------]\n")
    printf(" 1 -------- Methode de Jacobi -------------------\n")
    printf(" 2 ---------Methode de Gauss-Seidel -------------\n")
    printf(" 3 - -------Methode de Relaxation ---------------\n")
    choix=input("Choix : ")
    while choix<>1 & choix<>2 & choix<>3
          printf("Erreur de choix \n");
          printf(" 1 - Methode de Jacobi \n");
          printf(" 2 - Methode de Gauss-Seidel \n");
          printf(" 3 - Methode de Relaxation \n");
          choix=input("Choix : ")
    end


   switch choix
    case 1
        printf("METHODE DE JACOBI \n")
        // saisie de la taille de la matrice
        n= 25
        while n<=1 //verification
            printf("Erreur de taille ")
            n= input("Entrez une taille correct de la matrice : ")
        end
        // saisie des valeurs de la matrice
        diago=4;
        // verification pour avoir un determinant non nul
        while diago==0
            printf("Determinant de la matrice est nul \n")
            diago=input("Entrer l element de la diagonale non nul = ")
        end
        alpha=1;
        g_amma=1;
        [A,b,c]= tridiagAb(25,1,4,1)
        //saisi du terme d'erreur
        epsi=10^(-10)
        // saisie nombre d iterration max
        n_max_iter= n+1
        // verification pour avoir un d'iterration asser grand
        while n_max_iter <= n
            printf("Nombre d iterration maximum faible \n")
            n_max_iter= input("Nombre maximum d iterration : ")
        end
        [x,r_kj,erreurJ,rspec]=Jacobi(A,b,c,epsi,n_max_iter)
        // Tracé des vecteurs
        
    case 2
        printf("METHODE DE GAUSS-SEIDEL \n")
        // saisie de la taille de la matrice
        n= 25;
        while n<=1 //verification
            printf("Erreur de taille ")
            n= input("Entrez une taille correct de la matrice : ")
        end
        // saisie des valeurs de la matrice
        diago=4;
        // verification pour avoir un determinant non nul
        while diago==0
            printf("Determinant de la matrice est nul \n")
            diago=input("Entrer l element de la diagonale non nul = ")
        end
        alpha=1;
        g_amma=1;
        [A,b,c]= tridiagAb(25,1,4,1)
        //saisi du terme d'erreur
        epsi=10^(-10)
        // saisie nombre d iterration max
        n_max_iter= n+1
        // verification pour avoir un d'iterration asser grand
        while n_max_iter <= n
            printf("Nombre d iterration maximum faible \n")
            n_max_iter= input("Nombre maximum d iterration : ")
        end
        [x,r_k,erreurG,rspec]=GaussSeidel(A,b,c,10^(-10),n_max_iter)
        
        
    case 3
        printf("METHODE DE RELAXATION \n")
        // saisie de la taille de la matrice
        n= 25
        while n<=1 //verification
            printf("Erreur de taille ")
            n= input("Entrez une taille correct de la matrice : ")
        end
        // saisie des valeurs de la matrice
        diago=4;
        // verification pour avoir un determinant non nul
        while diago==0
            printf("Determinant de la matrice est nul \n")
            diago=input("Entrer l element de la diagonale non nul = ")
        end
        alpha=1;
        g_amma=1;
        [A,b,c]= tridiagAb(25,1,4,1)
        //saisie du terme d'erreur
        epsi=10^(-10)
        // saisie du parametre omega
        omega=1.5
        while omega==0
            printf("Erreur le parametre Omega doit être non nul \n")
            omega=input("Entrer le parametre Omega = ")
        end
        // saisie nombre d iterration max
        n_max_iter= n+1
        // verification pour avoir un d'iterration asser grand
        while n_max_iter <= n
            printf("Nombre d iterration maximum faible \n")
            n_max_iter= input("Nombre maximum d iterration : ")
        end
        [x,r_kR,erreurR,rspec]=Relaxation(A,b,c,epsi,omega,n_max_iter)
      
    end

  clf;
        plot([0,log(erreurJ)],"r");
        plot([0,log(erreurG)],"b");
        plot([0,log(erreurR)],"g");
        legend(["erreurs jacobi", "erreurs gauss","erreurs relaxation"],location="northwest");
        xlabel("Axe des k");
        ylabel("Axe des e(k)");
        title("SUITE DES ERREURS RELATIVES","fontsize",3)
        figure();
        plot([0,log(r_kj)],"r");
        plot([0,log(r_k)],"b");
        plot([0,log(r_kR)],"g");
        legend(["erreurs jacobi", "erreurs gauss","erreurs relaxation"],location="northwest");
        xlabel("Axe des k");
        ylabel("Axe des r(k)");
        title("SUITE DES RESIDUS RELATIFS","fontsize",3)
        figure();
    
    
    
    
    
    
    



