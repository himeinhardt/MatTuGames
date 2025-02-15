function [Sh_pGen] = airportShapley(standAlone,out_uc)
    %AIRPORTSHAPLEY Determines airport problem's Shapley value for a non-sorted vector
    %Input should be 2 vectors of te same size
    
    [standAlone_sort, ind] = sort(standAlone(:,1));
    n_Types = out_uc(ind,1);
    tot_Types = length(standAlone_sort);
    
    for k=1:tot_Types
        r(k,1)=sum(n_Types(k:tot_Types));
    end
    
    Sh_pGen_order = zeros(tot_Types,1);
    for i=1:tot_Types
        
        for k=1:i
            if k==1
                Sh_pGen_order(i,1) = Sh_pGen_order(i,1) + (standAlone_sort(k))/r(k,1)  ;
            else
                Sh_pGen_order(i,1) = Sh_pGen_order(i,1) + (standAlone_sort(k)-standAlone_sort(k-1))/r(k,1)   ;
            end
        end
    end
   
    unsorted = 1:length(standAlone(:,1));
    newInd(ind) = unsorted;
    
    Sh_pGen = Sh_pGen_order(newInd);
end

% %Sequential old
% num_GenOverSeg = zeros(num_outages,time_steps);
% aux = inf(num_outages,time_steps);
% uniqueSeg = zeros(num_outages,time_steps);
% for t=1:time_steps
%     for i=1:num_outages %considera "pedazos" ordenados
%         for j=1:num_outages
%             if AS_payment(j,t) >= AS_payment(i,t)
%                 num_GenOverSeg(i,t) = num_GenOverSeg(i,t) + out_uc(j,t); 
%                 %Sh_pGen(j,t) = Sh_pGen(j,t) + sort_AS_payment(i,t)%./sort_uc(i,t);
%             end
%             if AS_payment(j,t) < AS_payment(i,t)
%                 if AS_payment(i,t) - AS_payment(j,t) <= aux(i,t)
%                     uniqueSeg(i,t) = AS_payment(i,t) - AS_payment(j,t);
%                     aux(i,t) = uniqueSeg(i,t);
%                 end
%             end
%         end
%     end
% end
% Sh_pGen = zeros(num_outages,time_steps);
% Sh_pCluster = zeros(num_outages,time_steps);
% for t=1:time_steps
%     for i=1:num_outages %considera "pedazos" ordenados
%         for j=1:num_outages
%             if AS_payment(j,t) <= AS_payment(i,t)
%                 Sh_pGen(i,t) = Sh_pGen(i,t) + uniqueSeg(j,t)./num_GenOverSeg(j,t); 
%                 Sh_pCluster(i,t) = Sh_pGen(i,t)*out_uc(i,t);
%             end
%         end
%     end
% end