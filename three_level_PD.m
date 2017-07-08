tic;
warning off all
% three level phase diagram

pdelta1=4;  % detuning
pdelta2=-2.5;  % detuning
pomega=5;   % coupling
pg=0.5;    % coupling
% pmu=-0.4;  % chemical potential
pGamma1=0.01;  % atom decay rate
pGamma2=0.01;  % atom decay rate
pkappa=0.01;  % cavity decay rate
pz=4; % coordination number
% pk=0.3;  % hopping

Na = 9; % Fock basis truncation for photon
dima = Na+1; % dimension of photon Hilbert space
dims = 3; % dimension of atom Hilbert space
dimtot = 3*dima; % dimension of the total Hilbert space

Ia = speye(dima); % identity on photon subspace
Is = speye(3); % identity on atom subspace
Itot = speye(dimtot); % identity on the complete Hilbert space

a = spdiags(sqrt(0:Na)',1,dima,dima);
a = kron(a,Is);

v1 = Is(:,1); % |e> state of atom
v2 = Is(:,2); % |g1> state of atom
v3 = Is(:,3); % |g2> state of atom

s00 = kron(Ia,v1*v1'); % |e><e|
s11 = kron(Ia,v2*v2'); % |g1><g1|
s22 = kron(Ia,v3*v3'); % |g2><g2|
s10 = kron(Ia,v1*v2'); % |g1><e|
s20 = kron(Ia,v1*v3'); % |g2><e|
s12 = kron(Ia,v2*v3'); % |g1><g2|

%Build the term containing the detunings:
HDelta = pdelta1*s00+(pdelta1-pdelta2)*s22;
%the coupling terms:
Hcoupling = pg*(a'*s10+a*s10')+pomega*(s20+s20');
%from which we build up the atom Hamiltonian:
H3L = HDelta+Hcoupling;
% decay of atom channel Gamma_1:
LG1 = 0.5*pGamma1*(2*kron(conj(s10),s10)-kron(Itot,s10'*s10)-kron(s10.'*conj(s10),Itot));
% decay of atom channel Gamma_2:
LG2 = 0.5*pGamma2*(2*kron(conj(s20),s20)-kron(Itot,s20'*s20)-kron(s20.'*conj(s20),Itot));
% decay of photon channel kappa:
Lk = 0.5*pkappa*(2*kron(conj(a),a)-kron(Itot,a'*a)-kron(a.'*conj(a),Itot));

L3L = -1i*kron(Itot,H3L)+1i*kron(H3L.',Itot)+LG1+LG2+Lk; %total Liouvillian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain phase diagram by self-consistent calculation

% NC=0;
err=0.00001;  % accuracy
Lavg=11;  % largest average psi

% sta_lgk=0.0;slgk=-0.12;Nx=25+1;
% sta_pmu=0.1;smu=-0.032;Ny=35+1;
sta_lgk=log10(2);slgk=-0.0230103;Nx=100+1;
sta_pmu=0.1;smu=-0.00234;Ny=300+1;

dataPD=zeros(Nx*Ny,4);
dataOB=zeros(Nx*Ny,5);

for jpar=1:Ny
    
    jpar
    
    pmu=sta_pmu+smu*(jpar-1);
%     psi=5*(rand+1i*rand);  % initial order parameter
    
    for ipar=1:Nx
        
        ipar
        
        plgk=sta_lgk+slgk*(ipar-1);
        pk=10^plgk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% self-consistent calculation

        dpsi=1d30;
        psi=1.5*err;  % initial order parameter
        
        for i=1:Lavg

            if dpsi<=err
                
                dataPD((jpar-1)*Nx+ipar,:)=[i-1,pk,pmu,abs(psi)];  % data of phase diagram
                dataOB((jpar-1)*Nx+ipar,:)=real([i-1,pk,pmu,trace(a'*a*rhoS1),trace((a')*a*(a')*a*rhoS1)]); % evaluate <na>, <(Delta na)^2>, g2
                break
            
            else
                
                Navg=5*i-3;
                NC=0;
                series=ones(1,Navg)*psi;  % initial series of order parameter
                
                while dpsi>err && NC<=round(200+40*(i-1))
                    
                    NC=NC+1;

                    psi=abs(series(mod(NC-1,Navg)+1));

                    % mean field term:
                    Hmean = -pz*pk*(conj(psi)*a+psi*a');
                    % chemical potential:
                    Hchem = -pmu*(a'*a+s22+s00);
                    % total Liouvillian
                    L=L3L-1i*kron(Itot,Hchem+Hmean)+1i*kron((Hchem+Hmean).',Itot);

                    opts.tol = 1d-4; 
                    [rhoS1,lambda0] = eigs(L,1,'sm',opts); %find eigenvector with largest real part
                    % eigen0 = lambda0 %check that the eigenvalue is 0

                    rhoS1 = reshape(rhoS1,dimtot,dimtot); %reshape eigenvector into a matrix
                    rhoS1 = rhoS1/trace(rhoS1); %normalize

                    series(mod(NC,Navg)+1)=abs(trace(a*rhoS1));
                    psi=mean(series);
                    series(mod(NC,Navg)+1)=psi;

                    dpsi=abs(abs(series(mod(NC,Navg)+1))-abs(series(mod(NC-1,Navg)+1)));
                
                end
                
            end
            
        end
        
        if i==Lavg
            if dpsi<=err
                dataPD((jpar-1)*Nx+ipar,:)=[i,pk,pmu,abs(psi)];  % data of phase diagram
                dataOB((jpar-1)*Nx+ipar,:)=real([i,pk,pmu,trace(a'*a*rhoS1),trace((a')*a*(a')*a*rhoS1)]); % evaluate <na>, <(Delta na)^2>, g2
            else
                dataPD((jpar-1)*Nx+ipar,:)=[i,pk,pmu,100000];  % data of phase diagram
                dataOB((jpar-1)*Nx+ipar,:)=real([i,pk,pmu,100000,100000]); % evaluate <na>, <(Delta na)^2>, g2
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
end

save('dataPD.dat','dataPD','-ascii','-double')
save('dataOB.dat','dataOB','-ascii','-double')





% [V,D]=eig(full(L),'vector');
% rhoS2 = reshape(V(:,15),dimtot,dimtot); %reshape eigenvector into a matrix
% rhoS2 = rhoS2/trace(rhoS2); %normalize
% Pop2 = [trace(a'*a*rhoS2),trace(a*rhoS2)]; %evaluate populations
toc;