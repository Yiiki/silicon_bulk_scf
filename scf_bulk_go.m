%*************************************************************************%
% * Program : pseudopotential plane wave method with LDA                  %
% * Author : Yii                                                          %
% * Affiliation : Department of Microelectronics and Nanoelectronics, THU %
% * Date : 13-Jul-2020                                                    %
% ************************************************************************%
parpool('local',4)% launch parallel pool  delete(gcp('nocreate'))
%% 1 parameter definitions

%% (1) physical constants
har_eV = 27.2114;% Hartree , unit : eV
bohr = 5.2917725e-11;% Borh radius , unit: m

%% (2) geometry optimization
alta=(5.417:0.001:5.42).*1e-10;% conventional cell size, unit: \AA
alst=alta./bohr;% conventional cell size, unit: bohr
can=length(alst);% compute a0 number
rak_def=zeros(1,can);% rank deficient evaluation of cos(R.G) matrix

%% (3) crystal definitions
Z=4;% effective valence charge per ion , IV-group elements
blbv=[[.5,.5,0];[0,.5,.5];[.5,0,.5]]; % bravais lattice basis vector
tau_raw=1/8*[[-1 -1 -1];[1 1 1]]';% basis atoms definition
Fn=[2,2,2,2];% Energy occupation number
load('kpoints_set_means.mat')% kpst : fcc means kpoint set 
load('silicon_U_ps_s2p2.mat')% oprs : s,p,...,psEt, rc_list, r, drdy, yms
alpha_str=1.67085;% diamond Madelung constant ,from Martin
alpha_si=alpha_str*2^(1/3);% own-calculated , from J.C.Slate

%% (4) numerical algorithm definitions
scf_max = 15;% iteration steps = scf_max - 1
tepuc_mat=ones(scf_max-1,can);% total energy per unit cell matrix
mix=1;% charge-mixing scheme
E_cut = 30.5*0.5;% cut-off energy , unit: Hartree 
kpi = 2;% k-points numbers, kpi = 1,2,3,4 ,using 1, 32, 256, 2048 points; 

%% 2 radial non-local pseudopotential : $U_{ps}^{non-local}(r)$
ban=size(tau_raw,2);% basis atoms number
nE =length(Fn);% 4 valence band
kk =kpst{kpi,1};% specify k-points sample
nk=size(kk,1); % Number of k-points
k_wi=ones(nk,1)./nk;% wight of the k-points
nps=size(oprs,2)-5;% pseudopotential l = 0 ,1 , ... , num_ps -1
fpae=oprs(1,end-4);% free pseudo-atom energy
rc_list=oprs(1:nps,end-3);% rc cut index , for l = 0 , 1 , ...
rm=oprs(:,end-2);% radial mesh
drdy=oprs(1:end,end-1);% dr/dy on y-mesh
D=oprs(2,end)-oprs(1,end);% y-mesh acc
unlc_r=oprs(:,1:nps)+Z.*kron(1./rm,ones(1,nps));% non-local pp. U_ps + Z/r 
ups2zr=unlc_r;% residue non-local pp. out of rc , initial
for i_rc = 1:nps
    ups2zr(1:rc_list(i_rc),i_rc)=0;% rc cut out
end
unlc=unlc_r-ups2zr;% cut out non-local pp. that r>rc

%% 3 scf geometry optimization
fprintf('\t Completion: ');
showTimeToCompletion;startTime=tic;% showTimeToCompletion
load_p=parfor_progress(can);% parfor_progress
parfor i_alst = 1: can
%% 3.1 a0-dependent : $G_{mesh},R_{mesh}$
a0=alst(i_alst);% conventional cell size, unit: bohr
a=a0*blbv; % unit cell with specified lattice constant
Va=abs(det(a));% unit cell volume
tau=a0.*tau_raw; % atom positions, Descartes coordinate
R_a= (3*Va/(4*pi*Z*ban))^(1/3);% Wigner-Seitz radius of Si
gamma_Ewald=-alpha_si*Z^2/(2*R_a);% gamma_Ewald , unit: Ha
kmax=2*pi/a0;% first brillouin zone boundary
k=kmax.*kk; % first brillouin zone boundary with specified a0
b=2*pi*eye(3)/(a');% reciprocal lattice primitive vectors as rows vector
nGxyz=ceil(2*sqrt(2*E_cut./(b.^2*ones(3,1))));% E_cut-dependent nGmax
cx=(-nGxyz(1):nGxyz(1))';
cy=(-nGxyz(2):nGxyz(2))';
cz=(-nGxyz(3):nGxyz(3))';
P = [ kron(kron(cx,cy.^0),cz.^0) ...
      kron(kron(cx.^0,cy),cz.^0) ...
      kron(kron(cx.^0,cy.^0),cz) ];% reciprocal lattice points in 
                                        % primitive basis
G=P *b;% reciprocal lattice points in x ,y ,z unit vector
nsG = G.^2*ones(3,1);% norm square of G
G = G(nsG<2*E_cut,:);% Retain all reciprocal vectors within the cut-off
P = P(nsG<2*E_cut,:);% G point, irrelevant with a0
nG = size(G,1);% Number of G-vectors
% Generate the potential part of the Hamiltonian
GmG2=kron(ones(1,nG),G)-kron(ones(nG,1),reshape(G',1,3*nG));%G-G'(nGx3nG)
GmGt=kron(ones(1,nG),P)-kron(ones(nG,1),reshape(P',1,3*nG));%P-P'(nGx3nG)
% relatation between the list of (Gi-Gj) and G matrix
[uni_dG,~,ic_dG]=unique(reshape(reshape(GmGt',3*nG^2,1),3,nG^2)','rows');
ndG=size(uni_dG,1);%number of distinct Gi-Gj that G we sample can produce
rG=cell(1,ndG);
for i_ndG=1:ndG
    fic=find(ic_dG==i_ndG)-1;
    % rG gives the location of uni_dG in G matrix
    rG{i_ndG}=[(fic-mod(fic,nG))./nG+1,mod(fic,nG)+1];
end
min_rN=2*max(uni_dG)-1;% to avoid rank-deficient 
huG=ceil(0.5.*(min_rN + [2 2 2]));% 2*huG - 2 = N , N to be even
ax=union(linspace(-0.5,0,huG(1)),linspace(0,0.5,huG(1)))';% xmesh
ay=union(linspace(-0.5,0,huG(2)),linspace(0,0.5,huG(2)))';% ymesh
az=union(linspace(-0.5,0,huG(3)),linspace(0,0.5,huG(3)))';% zmesh
a_mesh = [ kron(kron(ax,ay.^0),az.^0) ...
           kron(kron(ax.^0,ay),az.^0) ...
           kron(kron(ax.^0,ay.^0),az)];% a-mesh
ndR=size(a_mesh,1);% a-mesh points amount
mid_R=0.5*(1+ndR);% mid-point of a_mesh
mid_G=0.5*(1+ndG);% mid_point of uni_dG
Rhaf=a_mesh(mid_R:end,:);% half of a_mesh
Ghaf=uni_dG(mid_G:end,:);% half of uni_dG
acs=cos(2.*pi.*Rhaf*Ghaf');% cos(R.G)
rak_def(i_alst)=rank(acs)-size(acs,2);% rank deficiency

%% 3.2 k-independent , charge-independent : $S,U_{ps}^{local}$
S=exp(-1i*(GmG2) *kron(eye(nG),tau))*...% structure factor
    kron(eye(nG),ones(ban,1));
SGList=sum(exp(-1i.*(uni_dG*b)*tau),2);% list version of structure factor
clG=4.*pi./((abs(GmG2).^2*kron(eye(nG),ones(3,1)))...
    +diag(Inf.*ones(nG,1)));% coulomb clG = 4 pi /norm($G_{i}$-$G_{j}$).^2
clGList=4.*pi./(((uni_dG*b).^2)*ones(3,1));% list version of clG
clGList(0.5*(1+ndG))=0;% elimilate singularity at Delta G = 0
UpsGG=-Z./Va.*clG;% Ups(G',G)

%% 3.3 k-dependent , charge-independent : $T,V_{ps}$
VG=zeros(nG,nG,nk); % Initialize potential part of the Hamiltionian
for i_k=1:nk
    Unc=zeros(nG,nG);% legendre poly. matrix
    kz=k(i_k,:);
    Gkz=G+kz;
    sqG=sum(Gkz.^2,2);
    noG=sqrt(sqG);
    gip=(Gkz*Gkz')./(noG*noG');
    gip(isnan(gip))=1;
    for i_l = 1:nps
        bsj=bj(i_l-1,rm*noG');
        bsj(isnan(bsj))=1*(i_l==1);
        lgd=lp(i_l-1,gip);
        Unc =Unc + 4.*pi.*(2.*(i_l-1)+1)./Va.*...
            (bsj'*((unlc(:,i_l).*rm.^2.*drdy.*D).*bsj)).*lgd;
    end
    U = S.*( UpsGG + Unc);% Ups for k_i
    T=0.5.*spdiags(sqG,0,nG,nG);% T for k_i
    VG(:,:,i_k)=T+U;
end

%% 3.4 scf charge iteration
% initializing
rhoG_mat=zeros(ndG,scf_max);% the rho_G list
rhoG_mat(0.5*(1+ndG),1)=ban*Z/Va;% free-electron initial guess
tepuc=gamma_Ewald.*ones(scf_max-1,1);% total energy per unit cell

% start scf 
for i=1:scf_max-1
    
E=zeros(nE,nk); % Initialize the eigenvalues
C=zeros(nG,nE,nk);% coefficients matrix

%% (1) $V_{coul}$
rhG=rhoG_mat(:,i);
hG=rhG(mid_G:end);
VcoulGG=clG;
for i_coul=1:ndG
    VcoulGG(sub2ind(size(VcoulGG),rG{i_coul}(:,1),rG{i_coul}(:,2)))=clG(...
        sub2ind(size(clG),rG{i_coul}(:,1),rG{i_coul}(:,2))).*rhG(i_coul);
end

%% (2) $V_{xc}$
hR=2.*cos(2.*pi.*(Rhaf*Ghaf'))*hG-hG(1);% half space charge density
VxcR=zeros(mid_R,1);
for i_vxcr=1:mid_R
    VxcR(i_vxcr)=Vc(hR(i_vxcr));
end
VxcR=VxcR+Vx(hR);
VxcG=acs\VxcR;
VxcGGL=[0.5.*flip(VxcG(2:end));VxcG(1);0.5.*VxcG(2:end)];
VxcGG=zeros(nG,nG);
for i_xc=1:ndG
VxcGG(sub2ind(size(VxcGG),rG{i_xc}(:,1),rG{i_xc}(:,2)))=VxcGGL(i_xc);
end

%% (3) $(T+V_{ps}+V_{coul}+V_{xc})\psi_{nk}=\epsilon_{nk}\psi_{nk}$
for i_k=1:nk
    [ef_M,ev_D]=eig(VG(:,:,i_k)+VcoulGG+VxcGG);
    ev0=real(diag(ev_D));
    [s_ev,id_ev]=sort(ev0);
    C(:,:,i_k)=ef_M(:,id_ev(1:nE));
    E(:,i_k)=real(s_ev(1:nE));
end
%% (4) $\rho^{i+1}(G)$
    cc=permute(C,[1 3 2]);
    rGk=zeros(ndG,size(cc,2));
for dG_idx=1:ndG
 for iE=1:length(Fn)
    rGk(dG_idx,:)=rGk(dG_idx,:)+(2*pi)^-3.*Fn(iE).*...
        sum(cc(rG{dG_idx}(:,1),:,iE).*conj(cc(rG{dG_idx}(:,2),:,iE)));
 end
end
rhG_new=(2*pi).^3/Va.*(rGk*k_wi);% rho(G)
rhoG_mat(:,i+1)=mix.*rhG_new + (1-mix).*rhoG_mat(:,i);% mixing scheme
%% (5) $E_{total}$
hGn=rhG_new(mid_G:end);% half rho(G)
hRn=2.*cos(2.*pi.*(Rhaf*Ghaf'))*hGn-hGn(1);% half charge density
ExcR=epc(hRn)+epx(hRn);% eps_xc(R)
excG=acs\ExcR;% eps_xc(G)
ep_GGL=[0.5.*flip(excG(2:end));excG(1);0.5.*excG(2:end)];% eps_xc list of G
tepuc(i)=tepuc(i)...
         + Fn*E*k_wi...                  % eigen energy value integral
         + Va.*(0.5.*rhG_new-rhG)'*(rhG_new.*clGList)...% Coulomb diff
         + Va.*rhG_new'*(ep_GGL - VxcGGL);            % Exc difference
end
tepuc_mat(:,i_alst)=real(tepuc);% E_total per unit cell
%---------5.3--------has------been------done------above-------------------%
load_p=parfor_progress;% update counting index
showTimeToCompletion(load_p/100,[],[],startTime);% update expectation
end

%% 4 total energy per unit cell convergence plot
[emin,emin_loc]=min(min(tepuc_mat));
figure % plot
subplot(1,2,1)
plot(alst.*bohr.*1e10,min(tepuc_mat),'LineWidth',1.75)
hold on
plot(alst(emin_loc).*bohr.*1e10,emin,'*','LineWidth',0.75)
hold off
xlabel('$a_{0}(\AA)$','interpreter','latex','FontSize',14)
ylabel('$E_{a}(Ha)$','interpreter','latex','FontSize',14)
text(.6,.8,['$E_{a,min}$ = ',num2str(emin),' Ha'],'Units',...
    'normalized','interpreter','latex','FontSize',14);
text(.6,.7,['with $a_{0}$ = ',num2str(alst(emin_loc).*bohr.*1e10),...
    ' $\AA$'],'Units','normalized','interpreter','latex','FontSize',14);
text(.6,.6,['and $E_{b}$ = ',num2str((fpae-emin*0.5).*har_eV),...
    ' $eV$'],'Units','normalized','interpreter','latex','FontSize',14);
subplot(1,2,2)
semilogy(abs(diff(tepuc_mat)),'LineWidth',1.75)
xlabel('$SCF$ $steps$','interpreter','latex','FontSize',14)
ylabel('$|E_{a}^{(i+1)}-E_{a}^{(i)}|$ $(Ha)$','interpreter','latex',...
    'FontSize',14)