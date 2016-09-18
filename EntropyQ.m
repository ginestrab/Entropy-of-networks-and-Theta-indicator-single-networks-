
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%  Entropy_Q (types) code %%%%%%%%%%%%%%%%%%%%%%
%This program evaluates the entropy of networks with given class assignment (type)
% The input of the program is a square adjacency matrix m of dimension n
% and a vector of dimension n describing the assignment type(i)=q
% size(type)=n,1
% with q integer between 1 and Q identifying Q classes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output of the program is the entropy Sigmac  
% Code part of the Supplementary material: Assessing the relevance of node features for network structure
% Ginestra Bianconi, Paolo Pin, and Matteo Marsili1, PNAS, April 2009
%
%
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This program is distributed ny the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%  
% If you use this code please cite 
%
% [1] G. Bianconi, P. Pin, M. Marsili
% Assessing the relevance of node features for network structure
% Proceedings of the National Academy of Sciences 106, no. 28 (2009): 11433-11438.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sigmac,W] = EntropyQ(m,type,n,Q)




precision = 10^(-3);
pi = 3.14159265358979323846;
undN = m;
undN = undN>0;
und = sum(undN)';
avgconn = sum(und)/n;

% Compute entropy with given degree sequence and community assignment --------
A=zeros(Q,Q);
for i=1:Q
   for j=1:Q
       A(i,j)=sum(sum( ( undN(type==i,type==j) ) )) ;
   end
end


% compute exp(Lagrangian multipliers)
z=rand(n,1); % z=exp(omega)
z=und/sum(und);
W=A/sum(sum(A));

% W(i,j)=exp( w(type(i),type(j)) )
oldW=zeros(Q,Q); oldz=zeros(n,1);
Scold=100;
for kk=1:2000
   bigW=W(type,type);
       U=(ones(n,1)*z' ).* bigW;
       D=ones(n,n) + ( z*z'.* bigW);
       U=max( U, 10^(-15) );
       D=max( D, 10^(-15) );
       z=und ./ (sum( ( U./D - diag(diag(U./D)) )' )' );
       z=max(z,10^(-15));
       z=z.*(und>0)+0*(und==0);
   for k=1:10
       d=zeros(Q,Q);
       bigW=W(type,type);
       for i=1:Q,
           for j=1:Q
               a=(type==i)+0; b=(type==j)+0;
              
               M=((a*b') .* (z*z') )./ ( ones(n,n) + (z*z' * W(i,j))) ;
              
          
               d(i,j)=( sum(sum( M-diag(diag(M)) )) );
                if(d(i,j)>0)
                   W(i,j)=A(i,j)/d(i,j);
                else
                    W(i,j)=0;
                end
           end
           
       end
   end
   
M=log( ones(n,n) + ( (z*z').*W(type,type) ) ) ;


z=z+0*(und==0);
Sc=(1/n)*( - sum(log(z+(z==0)).*und) - sum(sum( triu(A.* log(W+(W==0)) ) )) + sum(sum( triu( M,1) )) );
Sc=Sc-(1/n)*(0.5*sum(log(2*pi*und+(und==0)))-0.5*sum(sum(triu(log(2*pi*A+(A==0))))));

   if kk>150&&max(max(abs(Sc-Scold)))<precision 
       break
   end
   Scold=Sc;
end


M=log( ones(n,n) + ( (z*z').*W(type,type) ) ) ;


z=z+0*(und==0);
Sc=(1/n)*( - sum(log(z+(z==0)).*und) - sum(sum( triu(A.* log(W+(W==0)) ) )) + sum(sum( triu( M,1) )) );
Sc=Sc-(1/n)*(0.5*sum(log(2*pi*und+(und==0))))-0.5*sum(sum(triu(log(2*pi*A+(A==0))))));
Sigmac=Sc;



