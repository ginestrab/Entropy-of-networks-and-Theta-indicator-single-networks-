%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Entropy distance  code %%%%%%%%%%%%%%%%%%%%%%%%%% 
% The input of this function are the weighted undirected adjacency matrix pp
% The matrix dd of distances between each node i and node j of the network
% The number of bins Nd that we consider in the distance classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output are 
% I) the entropy Sigma_o of undirected networks with the
% same degree distribution of the matrix pp
% II) The entropy Sigma_d of the undirected networks with the same distribution
% of links at a given distance d of the network pp.
% III)The vector W(d) that modulates the probability that a node i is 
% connected with a node j at distance d from i 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sigma_0,Sigma_d,W]=Entropy_distance(pp,dd,Nd,l);
precision=10^(-2);
loops=1000;
pi=3.14159265358979323846;
Sigma_0=0;
Sigma_d=0;

    n=max(size(pp));
    Dmax=max(max(dd));
    [I,J,V]=find(dd);
    d_min=min(V);
   % l=logspace(log(d_min)/log(10),log(Dmax)/log(10),Nd);
   %l=linspace(-14,55,30);
 % l=linspace(0.2,9,Nd);
    for i=1:n,
        for j=1:n,
            h=hist(dd(i,j),l);
            I=find(h);
            class(i,j)=I(1);
        end
    end
    
	undN=+(pp>0);
    und=sum(undN)';
    avg_conn=sum(und)/n;
	
	% Compute no-partition entropy    ---------------------------------------------------------------------
	
	% compute exp(Lagrangean multipliers)
	z=rand(n,1);   % z=exp(omega)
    oldz=zeros(n,1);
	for kk=1:loops
        for k=1:10
            U=ones(n,1)*z';
            D=ones(n,n) + z*z';  
            z=und ./ (sum( ( U./D - diag(diag(U./D)) )' ) )';
            z=max(z,10^(-15));
        end
        if max(abs((z>0).*(1-z./(oldz+(oldz==0)))))<precision
            break
        end
        oldz=z;
	end
    M=log( ones(n,n)+z*z' );
    M2=(z*z')./( ( ones(n,n)+z*z' ).^2 );
    alpha=sum(M2);
	S=(1/n)*( - sum(log(z).*und) + sum(sum( triu(M,1) )) - ( sum(log(2*pi*alpha)) )/2 );
	Sigma_0=S;

    % Compute distance entropy    ---------------------------------------------------------------------
	B=zeros(Nd,1);
	for d=1:Nd
            B(d)=sum(sum(undN.*(class==d)))/2;
    end

	% compute exp(Lagrangean multipliers)
	z=rand(n,1);   % z=exp(omega)
	W=rand(Nd,1);    % W(d)=exp( w(d)) 
    oldW=zeros(Nd,1); oldz=zeros(n,1);
	for kk=1:loops
        bigW=zeros(n,n);
        for d=1:Nd,
        bigW=bigW+W(d)*(class==d);
        end
        U=ones(n,1)*z' .* bigW;
        D=ones(n,n) + ( z*z'.* bigW); 
        z=und ./ (sum( ( U./D - diag(diag(U./D)) )' ) )';
        z=max(z,10^(-15));
        
        B2=zeros(Nd,1);
        for d=1:Nd
                M=(class==d) .* (z*z') ./ ( ones(n,n) + ((z*z') .* bigW)) ;
                B2(d)=( sum(sum( (M)-diag(diag(M)) )) )/2;
                if (B2(d)*B(d))>0.0
                    W(d)=B(d)/(B2(d));
                    W(d)=max(W(d),10^(-15));
                    W(d)=min(W(d),10^15);
                else
                    W(d)=0;
                end
        end	
            
        if kk>30 && max(max(abs((W>0).*(1-W./(oldW+(oldW==0))))))<precision && max(abs((z>0).*(1-z./(oldz+(oldz==0)))))<precision
            break
        end
        oldW=W; oldz=z;
        end
     bigW=zeros(n,n);
        for d=1:Nd,
        bigW=bigW+W(d)*(class==d);
        end
    M=log( ones(n,n) + ( (z*z').*bigW ) ) ;
    
    Sc=(1/n)*( - sum(log(z).*und) - sum(B.*  log(W))  + sum(sum( triu( M,1) )) - ( sum(log(2*pi*und+(und==0))) )/2  - sum(log(2*pi*B+(B==0)))/2 );
    Sigma_d=Sc;
    



