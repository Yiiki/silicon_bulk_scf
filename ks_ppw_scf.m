function ksop=ks_ppw_scf(a0,Z,blbv,tau_raw,Fn,kpst,oprs,alpha_si,scf_max,mix,E_cut,kpi)
%*************************************************************************%
% * title : Kohn-Sham pseudopotential plane wave method with LDA          %
% * Author : Yii                                                          %
% * Affiliation : Department of Microelectronics and Nanoelectronics, THU %
% * Date : 18-Jul-2020                                                    %
%*************************************************************************%
ksop{7}=[];
if nargin==12
    % import post-process
    ban=size(tau_raw,2);% basis atoms number
    nE =length(Fn);% 4 valence band
    kk =kpst{kpi,1};% specify k-points sample
    nk=size(kk,1); % Number of k-points
    k_wi=ones(nk,1)./nk;% wight of the k-points
    nps=size(oprs,2)-4;% pseudopotential l = 0 ,1 , ... , num_ps -1
    unlc=oprs(:,1:nps);% rc cut pseudopotential l-th component
    pate=oprs(1,end-3);% pseudo atom total energy
    rm=oprs(:,end-2);% radial mesh
    drdy=oprs(:,end-1);% dr/dy on y-mesh
    D=oprs(2,end)-oprs(1,end);% y-mesh acc
    % a0-dependent 
    a=a0*blbv; % unit cell with specified lattice constant
    Va=abs(det(a));% unit cell volume
    tau=a0.*tau_raw; % atom positions, Descartes coordinate
    R_a= (3*Va/(4*pi*Z*ban))^(1/3);% Wigner-Seitz radius of Si
    gamma_Ewald=-alpha_si*Z^2/(2*R_a);% gamma_Ewald , unit: Ha
    kmax=2*pi/a0;% first brillouin zone boundary
    k=kmax.*kk; % first brillouin zone boundary with specified a0
    b=2*pi*eye(3)/(a');%reciprocal lattice primitive vectors as rows vector
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
    ndG=size(uni_dG,1);%number of distinct Gi-Gj
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

    % k-independent, charge-independent : $S,U_{ps}^{local}$
    S=exp(-1i*(GmG2) *kron(eye(nG),tau))*...% structure factor
        kron(eye(nG),ones(ban,1));
    clG=4.*pi./((abs(GmG2).^2*kron(eye(nG),ones(3,1)))...
        +diag(Inf.*ones(nG,1)));% coulomb clG = 4 pi /norm(Gi-Gj).^2
    clGList=4.*pi./(((uni_dG*b).^2)*ones(3,1));% list version of clG
    clGList(0.5*(1+ndG))=0;% elimilate singularity at Delta G = 0
    UpsGG=-Z./Va.*clG;% Ups(G',G)

    % scf charge iteration
    % initializing
    rhoG_mat=zeros(ndG,scf_max);% the rho_G list
    rhoG_mat(0.5*(1+ndG),1)=ban*Z/Va;% free-electron initial guess
    tepuc=gamma_Ewald.*ones(scf_max-1,1);% total energy per unit cell

    % start scf 
    for i=1:scf_max-1

    E=zeros(nE,nk); % Initialize the eigenvalues
    C=zeros(nG,nE,nk);% coefficients matrix

    % (1) $V_{coul}$
    rhG=rhoG_mat(:,i);
    hG=rhG(mid_G:end);
    VcoulGG=clG;
    for i_coul=1:ndG
     VcoulGG(sub2ind(size(VcoulGG),rG{i_coul}(:,1),rG{i_coul}(:,2)))...
     =clG(sub2ind(size(clG),rG{i_coul}(:,1),rG{i_coul}(:,2))).*rhG(i_coul);
    end

    % (2) $V_{xc}$
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

    % (3) $(T+V_{ps}+V_{coul}+V_{xc})\psi_{nk}=\epsilon_{nk}\psi_{nk}$
    for i_k=1:nk
        % 3.3 k-dependent , charge-independent : $T,V_{ps}$
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
        T=  0.5.*spdiags(sqG,0,nG,nG);% T for k_i
        % solve eig equation
        [ef_M,ev_D]=eig(T+U+VcoulGG+VxcGG);
        ev0=real(diag(ev_D));
        [s_ev,id_ev]=sort(ev0);
        C(:,:,i_k)=ef_M(:,id_ev(1:nE));
        E(:,i_k)=real(s_ev(1:nE));
    end
    % (4) $\rho^{i+1}(G)$
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
    % (5) $E_{total}$
    hGn=rhG_new(mid_G:end);% half rho(G)
    hRn=2.*cos(2.*pi.*(Rhaf*Ghaf'))*hGn-hGn(1);% half charge density
    ExcR=epc(hRn)+epx(hRn);% eps_xc(R)
    excG=acs\ExcR;% eps_xc(G)
    ep_GGL=[0.5.*flip(excG(2:end));excG(1);0.5.*excG(2:end)];% eps_xc list of G
    tepuc(i)=tepuc(i)...
             + Fn*E*k_wi...                  % eigen energy value integral
             + Va.*(0.5.*rhG_new-rhG)'*(rhG_new.*clGList)...% Coulomb diff
             + Va.*rhG_new'*(ep_GGL - VxcGGL)...          % Exc difference
             - ban*pate;                     % two pseudo free atom energy
    end
    ksop{1}=G;% G matrix
    ksop{2}=S;% structure factor
    ksop{3}=UpsGG;% pseudopotential
    ksop{4}=VcoulGG;% coulomb potential matrix
    ksop{5}=VxcGG;% xc energy matrix
    ksop{6}=rhG_new;% rho(G)
    ksop{7}=real(tepuc)./ban;% cohesive energy per atom

elseif nargin==0 % default parameters
    struct1=load('kpoints_set_means.mat');
    struct2=load('silicon_U_ps_s2p2.mat');
    har_eV = 27.2114;% Hartree , unit : eV
    bohr = 5.2917725e-11;% Borh radius , unit: m
    a0=5.42e-10/bohr;% 01 lattice constant
    Z=4;% 02 effective valence charge per ion , IV-group elements
    blbv=[[.5,.5,0];[0,.5,.5];[.5,0,.5]]; % 03 bravais lattice basis vector
    tau_raw=1/8*[[-1 -1 -1];[1 1 1]]';% 04 basis atoms definition
    Fn=[2,2,2,2];% 05 Energy occupation number
    kpst=struct1.kpst;% 06 kpst : fcc means kpoint set 
    oprs=struct2.oprs;% 07 oprs : s,p,...,psEt, rc_list, r, drdy, yms
    alpha_str=1.67085;% diamond Madelung constant ,from Martin
    alpha_si=alpha_str*2^(1/3);% 08 own-calculated , from J.C.Slate
    scf_max = 15;% 09 iteration steps = scf_max - 1
    mix=1;% 10 charge-mixing scheme
    E_cut = 30.5*0.5;% 11 cut-off energy , unit: Hartree 
    kpi = 2;%12 k-points numbers, kpi = 1,2,3,4 ,using 1, 32, 256, 2048 points; 
    ksop=ks_ppw_scf(a0,Z,blbv,tau_raw,Fn,kpst,oprs,alpha_si,scf_max,mix,E_cut,kpi);
    
    [emin,eloc]=min(ksop{end});
    figure % plot
    subplot(1,2,1)
    plot(1:(scf_max-1),ksop{end}.*har_eV,'-o','LineWidth',1.75)
    hold on
    plot(eloc,emin.*har_eV,'*','LineWidth',0.75)
    hold off
    xlabel('$SCF$ $Steps$','interpreter','latex','FontSize',14)
    ylabel('$-E_{b}(eV/atom)$','interpreter','latex','FontSize',14)
    text(.6,.8,['$E_{b}$ = ',num2str(-emin*har_eV),' eV'],'Units',...
    'normalized','interpreter','latex','FontSize',14);
    
    subplot(1,2,2)
    semilogy(abs(diff(ksop{end})),'LineWidth',1.75)
    xlabel('$SCF$ $steps$','interpreter','latex','FontSize',14)
    ylabel('$|E_{b}^{(i+1)}-E_{b}^{(i)}|$ $(Ha)$','interpreter','latex',...
    'FontSize',14)
end

end
