function [Nucl_pGen] = airportNucleolus(standAlone,out_uc)
%AIRPORTNUCLEOLUS Determines airport problem's Nucleolus for a non-sorted vector
%Input should be 2 vectors of te same size


    [standAlone_sort, ind] = sort(standAlone(:,1));
    n_Types = out_uc(ind,1);
    tot_Types = length(standAlone_sort);

    
    f = zeros(tot_Types,1); 
    alpha_q = zeros(tot_Types,1);                       %q depende del problema de optimizacion
    k_q = zeros(tot_Types,1);                           %q depende del problema de optimizacion
    min_q = zeros(tot_Types,1);                         %q depende del problema de optimizacion
    col_min_q = Inf(tot_Types, tot_Types);              %columna q depende del problema de optimizacion
    t_k = zeros(tot_Types, tot_Types);                  %columna q depende del problema de optimizacion
    
    q = 1;
    for k=1:tot_Types
        t_k(k,q) = sum(   n_Types( 1:k )   );           %se elimina k_q(q-1)+1 para este caso
    end
    
    for k=1:tot_Types-1
        col_min_q(k,q) = standAlone_sort(k,1)/(t_k(k,q)+1);
    end
    
    [min_q(q), k_q(q)] = min(col_min_q(:,q));            %hasta id_min f(1:id=min)=alpha 
    
    alpha_q(q) = -min(   min_q(q), standAlone_sort(tot_Types)/t_k(tot_Types,q)   )  ;
    alpha_q(q)
    
    if k_q(q) > tot_Types-2
        if k_q(q) == tot_Types
            f(1:k_q(q)) = -alpha_q(q);
        else   %tot_Types-1
            f(1:k_q(q)) = -alpha_q(q);
            f(tot_Types) = (standAlone_sort(tot_Types)+alpha_q(q)*t_k(tot_Types-1,q)) / n_Types(tot_Types);
        end
    else 
        f(1:k_q(q)) = -alpha_q(q);
        while  k_q(q) <= tot_Types-1 && alpha_q(q) == -min_q(q) 
    
            q=q+1;
            for k = k_q(q-1)+1 : tot_Types
                t_k(k,q) = sum(   n_Types( k_q(q-1)+1:k )   );
            end
    
            sum_aux = 0;
            for q_aux = 1:q-1
                sum_aux = sum_aux + t_k( k_q(q_aux), q_aux )  *  alpha_q(q_aux); 
            end
            
            for k = k_q(q-1)+1 : tot_Types-1
                col_min_q(k,q) = (   standAlone_sort(k,1) + sum_aux   )   /   ( t_k(k,q)+1 );  %falta completar suma para caso de q>=3
            end
            
            [min_q(q),k_q(q)] = min(col_min_q(:,q)); 
        
            alpha_q(q) = -min(   min_q(q),   (standAlone_sort(tot_Types,1)+ sum_aux)/t_k(tot_Types,q) );
            
            if alpha_q(q) == -min_q(q)
                f(k_q(q-1)+1:k_q(q)) = -alpha_q(q);    
            end
        end
    
        f(tot_Types) = (standAlone_sort(tot_Types,1)+ sum_aux)/t_k(tot_Types,q);
    end

    unsorted = 1:length(standAlone(:,1));
    newInd(ind) = unsorted;
    
    Nucl_pGen = f(newInd);
end

