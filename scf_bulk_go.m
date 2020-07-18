%*************************************************************************%
% * Program : pseudopotential plane wave method with LDA                  %
% * Author : Yii                                                          %
% * Affiliation : Department of Microelectronics and Nanoelectronics, THU %
% * Date : 13-Jul-2020                                                    %
% ************************************************************************%
parpool('local',4)% launch parallel pool  delete(gcp('nocreate'))

%% parameter definitions: (1) physical constants
har_eV = 27.2114;% Hartree , unit : eV
bohr = 5.2917725e-11;% Borh radius , unit: m

%%  parameter definitions: (2) geometry optimization
alta=(5.41:0.01:5.43).*1e-10;% conventional cell size, unit: \AA
alst=alta./bohr;% conventional cell size, unit: bohr
can=length(alst);% compute a0 number
%%  parameter definitions: (3) calculation definitions
Z=4;% effective valence charge per ion , IV-group elements
blbv=[[.5,.5,0];[0,.5,.5];[.5,0,.5]]; % bravais lattice basis vector
tau_raw=1/8*[[-1 -1 -1];[1 1 1]]';% basis atoms definition
Fn=[2,2,2,2];% Energy occupation number
load('kpoints_set_means.mat')% kpst : fcc means kpoint set 
load('silicon_U_ps_s2p2.mat')% oprs : s,p,...,psEt, r, drdy, yms
alpha_str=1.67085;% diamond Madelung constant ,from Martin
alpha_si=alpha_str*2^(1/3);% own-calculated , from J.C.Slate
scf_max = 15;% iteration steps = scf_max - 1
tepuc_mat=ones(scf_max-1,can);% total energy per unit cell matrix
mix=1;% charge-mixing scheme
E_cut = 30.5*0.5;% cut-off energy , unit: Hartree 
kpi = 2;% k-points numbers, kpi = 1,2,3,4 ,using 1, 32, 256, 2048 points; 

%% geometry optimization by brute force enumeration
fprintf('\t Completion: ');
% showTimeToCompletion;startTime=tic;% showTimeToCompletion
% load_p=parfor_progress(can);% parfor_progress
tic
parfor i_alst = 1: can % geometry optimization by brute force enumeration
    disp(i_alst);
a0=alst(i_alst);% conventional cell size, unit: bohr
ksop=ks_ppw_solver(a0,Z,blbv,tau_raw,Fn,kpst,oprs,alpha_si,scf_max,mix,...
    E_cut,kpi);% ksop={G, S, UpsGG, VcoulGG, VxcGG, rhG_new, tepuc}
tepuc_mat(:,i_alst)=ksop{end};% E_total per unit cell
%---------5.3--------has------been------done------above-------------------%
% load_p=parfor_progress;% update counting index
% showTimeToCompletion(load_p/100,[],[],startTime);% update expectation
end
toc

% output config
[emin,emin_loc]=min(min(tepuc_mat));
ap=alst(emin_loc);
ksbt=ks_ppw_solver(ap,Z,blbv,tau_raw,Fn,kpst,oprs,alpha_si,scf_max,mix,...
    E_cut,kpi);% ksbt={G, S, UpsGG, VcoulGG, VxcGG, rhG_new, tepuc}
gomi{1}=ap;% a0 unit: bohr
gomi{2}=ksbt{1};% G matrix
gomi{3}=ksbt{2};% structure factor
gomi{4}=ksbt{3};% pseudopotential
gomi{5}=ksbt{4};% coulomb potential matrix
gomi{6}=ksbt{5};% xc energy matrix
save(['Ecut_',num2str(2*E_cut),'Rdy_aG.mat'],'gomi')

% plot config
figure 
subplot(1,2,1)
plot(alst.*bohr.*1e10,min(tepuc_mat).*har_eV,'LineWidth',1.75)
hold on
plot(alst(emin_loc).*bohr.*1e10,emin.*har_eV,'*','LineWidth',0.75)
hold off
xlabel('$a_{0}(\AA)$','interpreter','latex','FontSize',14)
ylabel('$-E_{b}(eV)$','interpreter','latex','FontSize',14)
text(.6,.7,['$E_{b}$ = ',num2str(-emin.*har_eV),...
    ' $eV$'],'Units','normalized','interpreter','latex','FontSize',14);
text(.6,.6,['at $a_{0}$ = ',num2str(alst(emin_loc).*bohr.*1e10),...
    ' $\AA$'],'Units','normalized','interpreter','latex','FontSize',14);
subplot(1,2,2)
semilogy(abs(diff(tepuc_mat)),'LineWidth',1.75)
xlabel('$SCF$ $steps$','interpreter','latex','FontSize',14)
ylabel('$|E_{a}^{(i+1)}-E_{a}^{(i)}|$ $(Ha)$','interpreter','latex',...
    'FontSize',14)
