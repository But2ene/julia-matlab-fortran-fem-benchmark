function [ u, timings ] = fem_poisson( nx, ny )
% Solves uxx + uyy = 1 on a unit square with Q1 finite elements.
% Returns the computed solution vector in u, and the timimngs vector where:
%
% t_grid   = timings[1] time for grid generation
% t_ptr    = timings[2] time to calculate matrix pointers
% t_asm    = timings[3] time for system matrix assembly
% t_rhs    = timings[4] time for right hand side assembly
% t_bdr    = timings[5] time to set boundary conditions
% t_sparse = timings[6] time to convert to sparse matrix format
% t_sol    = timings[7] time for solution

% Copyright 2013-2016 Precise Simulation, Ltd.
% J.S. Hysing 161008.


tic
i_cub  = 2;
sfun   = 'sf_quad_Q1';
coef_a = [1 1];
coef_f = 1;

if( nargin<1 )
  nx = 16;
end
if( nargin<2 )
  ny = nx;
end
grid = rectgrid(nx,ny);
p    = grid.p;
c    = grid.c;
a    = [];
b    = grid.b;
t_grid = toc;
clear grid


tic
[vRowInds,vColInds,vAvals,n_rows,n_cols,t_ptr] = ...
          assemblea( [2 3;2 3], {sfun;sfun}, coef_a, i_cub, p, c, a );
t_asm = toc - t_ptr;


tic
ind_c  = b(1,:);
i_edge = b(2,:);
j_edge = mod( i_edge, size(c,1) ) + 1;
ix = sub2ind( size(c), [i_edge j_edge], [ind_c ind_c] );
indrow = unique( c(ix) );

pp = sparse( max(vRowInds), 1 ) ;
pp(indrow) = 1;
itmp = 1:length(vRowInds);
indr = itmp(logical(pp(vRowInds)));  % Index to entries in matrix corresponding to boundary rows.

vAvals(indr) = 0;   % Zero out bc rows.
indr = indr(1:length(indrow));
vRowInds(indr) = indrow;
vColInds(indr) = indrow;
vAvals(indr) = 1;   % Set diagonals to 1.
t_bdr = toc;
clear b


tic
A = sparse( vRowInds, vColInds, vAvals, n_rows, n_cols );
t_sparse = toc;
clear vRowInds vColInds vAvals n_rows n_cols


tic
[f,t_ptr_f] = assemblef( 1, {sfun}, coef_f, i_cub, p, c, a );
t_rhs = toc;
clear p c a


tic
f(indrow) = 0;
t_bdr = t_bdr + toc;


tic
u = A\f;
% u = zeros(size(f));
t_sol = toc;


timings = [ t_grid, t_ptr, t_asm, t_rhs, t_bdr, t_sparse, t_sol ];

function [ grid ]=rectgrid(gQFlFj,gQFlFjd,RJ7RK),if (~(nargin||nargout)),help rectgrid,return,end,if (nargin<3),RJ7RK=[0 1;0 1];end,if (nargin==1),gQFlFjd=gQFlFj;end,if (nargin==0),gQFlFj=10;gQFlFjd=10;end,RJ7RKZ=[];RJ7RKZ4=[];if (numel(gQFlFj)>1);RJ7RK(1)=min(gQFlFj);RJ7RK(3)=max(gQFlFj);RJ7RKZ=sort(gQFlFj);gQFlFj=numel(gQFlFj)-1;end,if (numel(gQFlFjd)>1);RJ7RK(2)=min(gQFlFjd);RJ7RK(4)=max(gQFlFjd);RJ7RKZ4=sort(gQFlFjd);gQFlFjd=numel(gQFlFjd)-1;end,Za3M_=RJ7RK(3)-RJ7RK(1);Za3M_H=RJ7RK(4)-RJ7RK(2);Za3M_HE=RJ7RK(1);PnF5K=RJ7RK(2);PnF5K2=gQFlFj*gQFlFjd;PnF5K2U=gQFlFj+1;EpPfN=gQFlFjd+1;EpPfNP=PnF5K2U*EpPfN;if (isempty(RJ7RKZ)),RJ7RKZ=linspace(Za3M_HE,Za3M_HE+Za3M_,PnF5K2U);end,if (isempty(RJ7RKZ4)),RJ7RKZ4=linspace(PnF5K,PnF5K+Za3M_H,EpPfN);end,[RJ7RKZ4,RJ7RKZ]=meshgrid(RJ7RKZ4,RJ7RKZ);grid.p=[ RJ7RKZ(:)'; RJ7RKZ4(:)' ];clear RJ7RKZ,clear RJ7RKZ4,EpPfNPT=reshape(1:EpPfNP,PnF5K2U,EpPfN);EpPfNPT=reshape(EpPfNPT(1:gQFlFj,1:gQFlFjd),1,PnF5K2);grid.c=[ EpPfNPT            ;EpPfNPT        + 1 ;EpPfNPT + PnF5K2U + 1 ;EpPfNPT + PnF5K2U     ];clear EpPfNPT,IUGvj=zeros(4,PnF5K2);IUGvj(1,:)=[ zeros(1,gQFlFj) 1:gQFlFj*(gQFlFjd-1) ];IUGvj(2,:)=[1:PnF5K2]+1;IUGvj(2,gQFlFj:gQFlFj:PnF5K2)=0;IUGvj(3,:)=[ gQFlFj+1:PnF5K2 zeros(1,gQFlFj) ];IUGvj(4,:)=[1:PnF5K2]-1;IUGvj(4,1:gQFlFj:PnF5K2-gQFlFj+1)=0;grid.a=IUGvj;clear IUGvj,IUGvjK(5,2*gQFlFj+2*gQFlFjd)=0;IUGvjK(1,:)=[  1:gQFlFj        gQFlFj:gQFlFj:gQFlFj*gQFlFjd  PnF5K2:-1:PnF5K2-gQFlFj+1  PnF5K2-gQFlFj+1:-gQFlFj:1 ];IUGvjK(2,:)=[  ones(1,gQFlFj)  2*ones(1,gQFlFjd)       3*ones(1,gQFlFj)     4*ones(1,gQFlFjd)     ];IUGvjK(3,:)=[  ones(1,gQFlFj)  2*ones(1,gQFlFjd)       3*ones(1,gQFlFj)     4*ones(1,gQFlFjd)     ];IUGvjK(4,:)=[ zeros(1,gQFlFj)    ones(1,gQFlFjd)        zeros(1,gQFlFj)      -ones(1,gQFlFjd)     ];IUGvjK(5,:)=[ -ones(1,gQFlFj)   zeros(1,gQFlFjd)         ones(1,gQFlFj)      zeros(1,gQFlFjd)     ];grid.b=IUGvjK;grid.s=ones(1,PnF5K2);
function [ T0VK7, T0VK7K, T0VK7Kb, NNX70l ]=sf_quad_Q1(g3E2m,g3E2mI,g3E2mIN,bMxdE,bMxdEl,bMxdElh,T0VK7),T0VK7K=[4 0 0 0];T0VK7Kb=[-1  1 1 -1;-1 -1 1  1];NNX70l='sf_quad_Q1';switch g3E2m,case 1,switch bMxdE,case 1,T0VK7=(1-bMxdEl(1))*(1-bMxdEl(2))/4;case 2,T0VK7=(1+bMxdEl(1))*(1-bMxdEl(2))/4;case 3,T0VK7=(1+bMxdEl(1))*(1+bMxdEl(2))/4;case 4,T0VK7=(1-bMxdEl(1))*(1+bMxdEl(2))/4;end,case {2,3},switch bMxdE,case 1,kilcN=-(1-bMxdEl(2))/4;kilcNr=-(1-bMxdEl(1))/4;case 2,kilcN=(1-bMxdEl(2))/4;kilcNr=-(1+bMxdEl(1))/4;case 3,kilcN=(1+bMxdEl(2))/4;kilcNr=(1+bMxdEl(1))/4;case 4,kilcN=-(1+bMxdEl(2))/4;kilcNr=(1-bMxdEl(1))/4;end,if (g3E2m==2),T0VK7=bMxdElh(:,1)*kilcN+bMxdElh(:,2)*kilcNr;elseif (g3E2m==3),T0VK7=bMxdElh(:,3)*kilcN+bMxdElh(:,4)*kilcNr;end,otherwise,T0VK7=0;end,function [ cVoBQ, cVoBQv, cVoBQvT, tQRBb, tQRBbu, tQRBbuA ]=assemblea(kilcNrH,JFHCh,JFHChE,NNX70,JfHH3E,JfHH3Ez,IUGvj),sind=ones(1,size(JfHH3Ez,2));JFHChE3=25000;SqmJO=[0 0];SqmJOA=[];g3E2mI=size(JfHH3E,1);g3E2mIN=size(JfHH3Ez,1);SqmJOAm=size(kilcNrH,2);Xfgik=g3E2mI*g3E2mIN;XfgikM=SqmJO(1);if (ischar(JFHCh)),JFHCh={JFHCh; JFHCh};elseif (numel(JFHCh)==1),JFHCh=[JFHCh; JFHCh];end,[~,XfgikMJ,~,c7TUd]=evalsfun(JFHCh{1},0,g3E2mI,g3E2mIN);[~,c7TUd3,~,c7TUd37]=evalsfun(JFHCh{2},0,g3E2mI,g3E2mIN);XfgikMJ=sum(XfgikMJ(:));c7TUd3=sum(c7TUd3(:));if (ischar(JFHChE)),JFHChE={ JFHChE };elseif (~iscell(JFHChE)),JFHChE=mat2cell(JFHChE,ones(1,size(JFHChE,1)),ones(1,size(JFHChE,2)));end,Z3XGk=numel(JFHCh)*2;Z3XGkz=false;Z3XGkzi=asmgroups(JFHChE3,sind);b8ssu=[];b8ssuD=[];tic,[cVoBQ,cVoBQv,tQRBb,tQRBbu]=asmrowcolptr(JFHCh,JfHH3E,JfHH3Ez,IUGvj,Z3XGkzi(2,:),b8ssu,b8ssuD,SqmJO);tQRBbuA=toc;b8ssuDO=4;b2J3v=0.577350269189626;bMxdEl=[ -b2J3v -b2J3v; b2J3v -b2J3v; -b2J3v b2J3v; b2J3v b2J3v ]';b2J3vC=ones(1,4);cVoBQvT=cell(XfgikMJ*c7TUd3,size(Z3XGkzi,2));for b2J3vCd=1:size(Z3XGkzi,2),iPu5h=length(Z3XGkzi{2,b2J3vCd});for iPu5h_=1:XfgikMJ*c7TUd3,cVoBQvT{iPu5h_,b2J3vCd}=zeros(iPu5h,1);end,end,iPu5h_H=[];K_rAD=[];K_rADV=[];K_rADVZ=[];fkvdG=0;for b2J3vCd=1:size(Z3XGkzi,2),fkvdGu=Z3XGkzi{1,b2J3vCd};fkvdGuk=Z3XGkzi{2,b2J3vCd};iPu5h=length(fkvdGuk);if (size(iPu5h_H,1)~=iPu5h),iPu5h_H=zeros(iPu5h,XfgikMJ);K_rAD=zeros(iPu5h,c7TUd3);end,[K_rADV,~]=tfjacquad(0,JfHH3E,JfHH3Ez(:,fkvdGuk),K_rADV,[],[],Z3XGkz);Xfgik=XfgikM;for k6QD6=1:b8ssuDO,[K_rADV,K_rADVZ]=tfjacquad(2,JfHH3E,JfHH3Ez(:,fkvdGuk),K_rADV,bMxdEl(:,k6QD6),K_rADVZ);if (k6QD6==1&&any(K_rADVZ(:,end)<=0)),error('assemblea: cell orientation wrong, negative determinant of Jacobian.'),end,for k6QD6g=1:SqmJOAm,k6QD6gU=JFHChE{min(k6QD6g,length(JFHChE))};if (ischar(k6QD6gU)),k6QD6gU=LJlLs(k6QD6gU,bMxdEl(:,k6QD6),fkvdGu,fkvdGuk,[],SqmJOA,K_rADVZ);end,LJlLsi=k6QD6gU.*K_rADVZ(:,end)*b2J3vC(k6QD6);for GIOIS=1:c7TUd3,K_rAD(:,GIOIS)=feval(c7TUd37,kilcNrH(2,k6QD6g),g3E2mI,g3E2mIN,GIOIS,bMxdEl(:,k6QD6),K_rADVZ,K_rAD(:,GIOIS));end,if (strcmp(c7TUd,c7TUd37)&&kilcNrH(1,k6QD6g)==kilcNrH(2,k6QD6g)),iPu5h_H=K_rAD;else,for GIOISG=1:XfgikMJ,iPu5h_H(:,GIOISG)=feval(c7TUd,kilcNrH(1,k6QD6g),g3E2mI,g3E2mIN,GIOISG,bMxdEl(:,k6QD6),K_rADVZ,iPu5h_H(:,GIOISG));end,end,for GIOIS=1:c7TUd3,for GIOISG=1:XfgikMJ,iPu5h_=(GIOIS-1)*XfgikMJ+GIOISG;cVoBQvT{iPu5h_,b2J3vCd}=cVoBQvT{iPu5h_,b2J3vCd}+iPu5h_H(:,GIOISG).*K_rAD(:,GIOIS).*LJlLsi;end,end,end,XfgikM=XfgikM+Z3XGk;end,XfgikM=XfgikM+1;end,cVoBQvT=cat(1,cVoBQvT{:});
function [ M1YtH, tQRBbuA ]=assemblef(kilcNrH,JFHCh,JFHChE,NNX70,JfHH3E,JfHH3Ez,IUGvj),sind=ones(1,size(JfHH3Ez,2));JFHChE3=25000;SqmJOA=[];g3E2mI=size(JfHH3E,1);g3E2mIN=size(JfHH3Ez,1);SqmJOAm=length(kilcNrH);Xfgik=g3E2mI*g3E2mIN;XfgikM=JFHChE3(1);if (iscell(JFHCh)),JFHCh=JFHCh{1};end,[~,LJlLsim,~,NNX70l]=evalsfun(JFHCh,0,g3E2mI,g3E2mIN);Z3XGk=numel(JFHCh)*XfgikM;Z3XGkz=false;tic,[tQRBb,TeodB]=mapdofbdr(JfHH3E,JfHH3Ez,IUGvj,[],LJlLsim);TeodB=TeodB';LJlLsim=sum(LJlLsim(:));tQRBbuA=toc;M1YtH=zeros(tQRBb,1);if (~iscell(JFHChE)),if (ischar(JFHChE)),JFHChE={ JFHChE };else,JFHChE=mat2cell(JFHChE,ones(1,size(JFHChE,1)),ones(1,size(JFHChE,2)));end,end,TeodBa=0;for GIOIS=1:length(JFHChE),if (~((ischar(JFHChE{GIOIS})&&strcmp(JFHChE{GIOIS},'0'))||(isnumeric(JFHChE{GIOIS})&&(JFHChE{GIOIS}==0)))),TeodBa=TeodBa+1;end,end,if (~TeodBa||isempty(kilcNrH)||~any(kilcNrH)),return,end,b8ssuDO=4;b2J3v=0.577350269189626;bMxdEl=[ -b2J3v -b2J3v; b2J3v -b2J3v; -b2J3v b2J3v; b2J3v b2J3v ]';b2J3vC=ones(1,4);T0VK7=[];K_rADV=[];K_rADVZ=[];Z3XGkzi=asmgroups(JFHChE3,sind);for b2J3vCd=1:size(Z3XGkzi,2),fkvdGu=Z3XGkzi{1,b2J3vCd};fkvdGuk=Z3XGkzi{2,b2J3vCd};iPu5h=length(fkvdGuk);[K_rADV,~]=tfjacquad(0,JfHH3E,JfHH3Ez(:,fkvdGuk),K_rADV,[],[],Z3XGkz);Xfgik=XfgikM;for k6QD6=1:b8ssuDO,[K_rADV,K_rADVZ]=tfjacquad(2,JfHH3E,JfHH3Ez(:,fkvdGuk),K_rADV,bMxdEl(:,k6QD6),K_rADVZ);if (k6QD6==1&&any(K_rADVZ(:,end)<=0)),error('assemblef: cell orientation wrong, negative determinant of Jacobian.'),end,for k6QD6g=1:SqmJOAm,k6QD6gU=JFHChE{min(k6QD6g,length(JFHChE))};if (ischar(k6QD6gU)),k6QD6gU=LJlLs(k6QD6gU,bMxdEl(:,k6QD6),fkvdGu,fkvdGuk,[],SqmJOA);end,TeodBa5=k6QD6gU.*K_rADVZ(:,end)*b2J3vC(k6QD6);for GIOIS=1:LJlLsim,T0VK7=evalsfun(NNX70l,kilcNrH(k6QD6g),g3E2mI,g3E2mIN,GIOIS,bMxdEl(:,k6QD6),K_rADVZ,T0VK7);jfsCg=[TeodB(fkvdGuk,GIOIS); tQRBb];M1YtH=M1YtH+accumarray(jfsCg,[T0VK7.*TeodBa5; 0]);end,end,XfgikM=XfgikM+Z3XGk;end,XfgikM=XfgikM+1;end,Z3XGk=Xfgik+XfgikM;
function [ K_rADV, K_rADVZ ]=tfjacquad(jfsCg_,JfHH3E,JfHH3Ez,K_rADV,bMxdEl,K_rADVZ,Z3XGkz),persistent jfsCg_u,if (nargin<7),Z3XGkz=false;end,if (jfsCg_==0||jfsCg_==1),PnF5K2=size(JfHH3Ez,2);if (~isequal(size(K_rADV),[PnF5K2 8])),K_rADV=zeros(PnF5K2,8);end,odscF=reshape(JfHH3E(1,JfHH3Ez'),PnF5K2,4);odscFV=reshape(JfHH3E(2,JfHH3Ez'),PnF5K2,4);K_rADV(:,1)=odscFV(:,3)-odscFV(:,2);K_rADV(:,2)=odscFV(:,4)-odscFV(:,1);K_rADV(:,3)=odscFV(:,1)-odscFV(:,2);K_rADV(:,4)=odscFV(:,4)-odscFV(:,3);K_rADV(:,5)=odscF(:,2)-odscF(:,3);K_rADV(:,6)=odscF(:,1)-odscF(:,4);K_rADV(:,7)=odscF(:,2)-odscF(:,1);K_rADV(:,8)=odscF(:,3)-odscF(:,4);if (Z3XGkz),jfsCg_u=K_rADV;end,end,if (jfsCg_>=1),PnF5K2=size(JfHH3Ez,2);if (~isequal(size(K_rADVZ),[PnF5K2 5])),K_rADVZ=zeros(PnF5K2,5);end,odscFVN=(1-bMxdEl(1))/4;eeX6H=(1+bMxdEl(1))/4;eeX6Hg=(1-bMxdEl(2))/4;eeX6Hg1=(1+bMxdEl(2))/4;K_rADVZ(:,1)=eeX6H*K_rADV(:,1)+odscFVN*K_rADV(:,2);K_rADVZ(:,2)=eeX6Hg*K_rADV(:,3)+eeX6Hg1*K_rADV(:,4);K_rADVZ(:,3)=eeX6H*K_rADV(:,5)+odscFVN*K_rADV(:,6);K_rADVZ(:,4)=eeX6Hg*K_rADV(:,7)+eeX6Hg1*K_rADV(:,8);K_rADVZ(:,5)=K_rADVZ(:,1).*K_rADVZ(:,4)-K_rADVZ(:,2).*K_rADVZ(:,3);K_rADVZ(:,1)=K_rADVZ(:,1)./K_rADVZ(:,5);K_rADVZ(:,2)=K_rADVZ(:,2)./K_rADVZ(:,5);K_rADVZ(:,3)=K_rADVZ(:,3)./K_rADVZ(:,5);K_rADVZ(:,4)=K_rADVZ(:,4)./K_rADVZ(:,5);end,if (~isempty(jfsCg_u)),K_rADV=jfsCg_u;end,function [ qrUac, TeodB, qrUacP ]=mapdofbdr(varargin),bMxdEl=[];if (nargin==5),qrUacPL=false;JfHH3E=varargin{1};JfHH3Ez=varargin{2};IUGvj=varargin{3};IUGvjK=varargin{4};T0VK7K=varargin{5};g3E2mI=size(JfHH3E,1);EpPfNP=size(JfHH3E,2);BTSXr=size(JfHH3Ez,1);elseif (nargin==3||nargin==4),qrUacPL=true;JfHH3Ez=varargin{1};IUGvjK=varargin{2};T0VK7K=varargin{3};if (nargin==4),bMxdEl=varargin{4};end,g3E2mI=2;BTSXr=size(JfHH3Ez,1);EpPfNP=max(JfHH3Ez(:));end,PnF5K2=size(JfHH3Ez,2);BTSXrR=(g3E2mI==2)*BTSXr+(g3E2mI==3)*(6*(BTSXr==4)+12*(BTSXr==8));BTSXrRI=(g3E2mI==3)*(BTSXr-2*(BTSXr==8));if (any(T0VK7K(:))<=0||any(rem(T0VK7K(:),1))||any(mod(T0VK7K(:,1),BTSXr)~=0)||any(mod(T0VK7K(:,2),BTSXrR)~=0)||any(mod(T0VK7K(:,3),BTSXrRI)~=0)),error('nLDof specification incorrect.'),end,tdXuM=sum(T0VK7K(:,1));tdXuMU=sum(T0VK7K(:,2));tdXuMUg=sum(T0VK7K(:,3));gIE5Y=sum(T0VK7K(:,4));gIE5YT=size(T0VK7K,1);TeodB=zeros(sum(T0VK7K(:)),PnF5K2);if (qrUacPL),qrUacP=cell(gIE5YT,1);[qrUacP{:}]=deal([]);else,qrUacP=[];end,gIE5YTd=nargout>2&&~isempty(IUGvjK);if (isempty(bMxdEl)),bMxdEl=[-1 1 1 -1;-1 -1 1 1];end,if (size(bMxdEl,2)~=sum(T0VK7K(:))),error('Wrong number of local cell coordinates.'),end,if (tdXuM>0),for dfY8Q=1:gIE5YT,if (T0VK7K(dfY8Q,1)>0),TeodB=update_dof_map(TeodB,JfHH3Ez);end,end,if (gIE5YTd),dfY8Qi=local_cell_dof_numbering(0,g3E2mI,BTSXr);dfY8QiT=compute_boundary_map(IUGvjK,dfY8Qi,bMxdEl,TeodB);for dfY8Q=1:gIE5YT,if (T0VK7K(dfY8Q,1)>0),if (qrUacPL),qrUacP{dfY8Q}=[ qrUacP{dfY8Q} dfY8QiT ];else,qrUacP=[ qrUacP dfY8QiT ];end,dfY8QiT(4,:)=dfY8QiT(4,:)+EpPfNP;end,end,end,end,if (tdXuMU>0),[fhOuO,fhOuOU]=fhOuOUK(JfHH3Ez,g3E2mI);MoT03=max([0; find(TeodB(:,1))]);for dfY8Q=1:gIE5YT,if (T0VK7K(dfY8Q,2)>0),MoT03f=T0VK7K(dfY8Q,2)/BTSXrR;MoT03fx=compute_edge_face_dof_map(1,g3E2mI,JfHH3Ez,fhOuO,fhOuOU,MoT03f);TeodB=update_dof_map(TeodB,MoT03fx);end,end,if (gIE5YTd),sQLYx=local_cell_dof_numbering(1,g3E2mI,BTSXr);sQLYxm=size(fhOuOU,1);for dfY8Q=1:gIE5YT,if (T0VK7K(dfY8Q,2)>0),MoT03f=T0VK7K(dfY8Q,2)/BTSXrR;sQLYxmf=[];for GIOISG=1:MoT03f,sQLYxmf=[ sQLYxmf sQLYx+MoT03 ];MoT03=MoT03+BTSXrR;end,NmyIn=compute_boundary_map(IUGvjK,sQLYxmf,bMxdEl,TeodB);if (qrUacPL),qrUacP{dfY8Q}=[ qrUacP{dfY8Q} NmyIn ];else,qrUacP=[ qrUacP NmyIn ];end,end,end,end,end,if (tdXuMUg>0),[M1YtH,NmyInO]=NmyInOE(JfHH3Ez);MoT03=max([0; find(TeodB(:,1))]);for dfY8Q=1:gIE5YT,if (T0VK7K(dfY8Q,3)>0),hXgro=T0VK7K(dfY8Q,3)/BTSXrRI;hXgroH=compute_edge_face_dof_map(2,g3E2mI,JfHH3Ez,M1YtH,NmyInO,hXgro);TeodB=update_dof_map(TeodB,hXgroH);end,end,if (gIE5YTd),hXgroHf=local_cell_dof_numbering(2,g3E2mI,BTSXr);uggXr=size(NmyInO,1);for dfY8Q=1:gIE5YT,if (T0VK7K(dfY8Q,3)>0),hXgro=T0VK7K(dfY8Q,3)/BTSXrRI;uggXrL=[];for GIOISG=1:hXgro,uggXrL=[ uggXrL hXgroHf+MoT03 ];MoT03=MoT03+BTSXrRI;end,uggXrL9=compute_boundary_map(IUGvjK,uggXrL,bMxdEl,TeodB);if (qrUacPL),qrUacP{dfY8Q}=[ qrUacP{dfY8Q} uggXrL9 ];else,qrUacP=[ qrUacP uggXrL9 ];end,end,end,end,end,if (gIE5Y>0),for dfY8Q=1:gIE5YT,for GIOISG=1:T0VK7K(dfY8Q,4),TeodB=update_dof_map(TeodB,[1:PnF5K2]);end,end,end,qrUac=max(TeodB(:));
function [ vRVBO ]=update_dof_map(vRVBO,vRVBOg),vRVBOg7=max([ 0; find( vRVBO(:,1) ) ]);qrUac=max(vRVBO(:));vRVBO(vRVBOg7+1:vRVBOg7+size(vRVBOg,1),:)=vRVBOg+qrUac;
function [ jqYEV ]=compute_edge_face_dof_map(jqYEVz,g3E2mI,JfHH3Ez,jqYEVzq,bQBCM,bQBCMd),if (bQBCMd==1),jqYEV=jqYEVzq;return,end,BTSXr=size(JfHH3Ez,1);PnF5K2=size(JfHH3Ez,2);if (jqYEVz==1&&g3E2mI==3),if (BTSXr==4),bQBCMdt=[ 1 2; 2 3; 3 1; 1 4; 2 4; 3 4 ];elseif (BTSXr==8),bQBCMdt=[ 1 2; 2 3; 3 4; 4 1; 1 5; 2 6; 3 7; 4 8; 5 6; 6 7; 7 8; 8 5 ];end,else,bQBCMdt=local_cell_dof_numbering(0,g3E2mI,BTSXr);end,NRhhf=size(bQBCMdt,2);NRhhfa=size(bQBCMdt,1);NRhhfaA=ones(size(jqYEVzq));for HTHgN=1:NRhhfa,HTHgNh=JfHH3Ez(bQBCMdt(HTHgN,:),:)';HTHgNhi=bQBCM(jqYEVzq(HTHgN,:),:);eiklt=find(any(HTHgNh-HTHgNhi,2));if (~isempty(eiklt)),for GIOIS=1:numel(eiklt),eikltt=eiklt(GIOIS);eiklttn=find(HTHgNh(eikltt,1)==HTHgNhi(eikltt,:));if (jqYEVz==1),eiklttn=-bQBCMd;elseif (HTHgNhi(eikltt,mod(eiklttn,size(bQBCMdt,2))+1)~=HTHgNh(eikltt,2)),eiklttn=-eiklttn;end,NRhhfaA(HTHgN,eikltt)=eiklttn;end,end,end,if (jqYEVz==2&&bQBCMd>NRhhf),m5RIl=find(mod(1:bQBCMd,floor(bQBCMd/NRhhf))==1);m5RIl=m5RIl(1:NRhhf)';NRhhfaA=sign(NRhhfaA).*m5RIl(abs(NRhhfaA));end,m5RIlm=false;if (jqYEVz==2&&g3E2mI==3),if (NRhhfa==6&&bQBCMd==9),m5RIlm=true;m5RIlme=zeros(18,9);m5RIlme([1 3 5 7 10 12 14 16],:)=[ 1 2 3 4 5 6 7 8 9 ;3 4 5 6 7 8 1 2 9 ;5 6 7 8 1 2 3 4 9 ;7 8 1 2 3 4 5 6 9 ;1 8 7 6 5 4 3 2 9 ;3 2 1 8 7 6 5 4 9 ;5 4 3 2 1 8 7 6 9 ;7 6 5 4 3 2 1 8 9 ];elseif (NRhhfa==6&&bQBCMd==16),m5RIlm=true;m5RIlme=zeros(32,16);m5RIlme([1 5 9 13 17 21 25 29],:)=[ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ;4 5 6 7 8 9 10 11 12 1 2 3 14 15 16 13 ;7 8 9 10 11 12 1 2 3 4 5 6 15 16 13 14 ;10 11 12 1 2 3 4 5 6 7 8 9 16 13 14 15 ;1 12 11 10 9 8 7 6 5 4 3 2 13 16 15 14 ;4 3 2 1 12 11 10 9 8 7 6 5 14 13 16 15 ;7 6 5 4 3 2 1 12 11 10 9 8 15 14 13 16 ;10 9 8 7 6 5 4 3 2 1 12 11 16 15 14 13 ];end,end,jqYEV(NRhhfa*bQBCMd,PnF5K2)=0;USJtk=max(jqYEVzq(:));for GIOIS=1:bQBCMd,if (m5RIlm),USJtkM=jqYEVzq+USJtk*(reshape(m5RIlme(abs(NRhhfaA)+bQBCMd*(NRhhfaA<0),GIOIS),size(NRhhfaA))-1);else,USJtkM=jqYEVzq+USJtk*(abs(NRhhfaA)-1);end,USJtkM_=(GIOIS-1)*NRhhfa+1;ipzeQ=GIOIS*NRhhfa;jqYEV(USJtkM_:ipzeQ,:)=USJtkM;if (~m5RIlm),NRhhfaA=next_values(NRhhfaA,bQBCMd);end,end,function [ IUGvjK ]=next_values(IUGvj,ipzeQz),ipzeQz0=sign(IUGvj)<0;IUGvjK=mod(IUGvj+ipzeQz0,ipzeQz)+(~ipzeQz0)-ipzeQz*ipzeQz0;
function [ uMjyj ]=compute_boundary_map(IUGvjK,uMjyjV,bMxdEl,TeodB),uMjyjV4=size(IUGvjK,2);rQrH1=size(uMjyjV,1);rQrH1L=size(uMjyjV,2);uMjyj(5+size(bMxdEl,1),uMjyjV4*rQrH1L)=0;rQrH1Lu=0;for oCAvJ=1:rQrH1,oCAvJU=find(IUGvjK(2,:)==oCAvJ);oCAvJUA=numel(oCAvJU);vhbrt=IUGvjK(1,oCAvJU);for GIOISG=1:rQrH1L,vhbrtO=uMjyjV(oCAvJ,GIOISG);uMjyj(:,rQrH1Lu+1:rQrH1Lu+oCAvJUA)=[ IUGvjK(1:3,oCAvJU);TeodB( vhbrtO, vhbrt );repmat( [      vhbrtO  ;bMxdEl(:,vhbrtO) ], 1, oCAvJUA ) ];rQrH1Lu=rQrH1Lu+oCAvJUA;end,end,function [ uMjyjV ]=local_cell_dof_numbering(vhbrtO9,varargin),switch (vhbrtO9),case 0,uMjyjV=[ 1 2 ;2 3 ;3 4 ;4 1 ];case 1,uMjyjV=[ 1:4 ]';end,function [ cVoBQ, cVoBQv, tQRBb, tQRBbu ]=asmrowcolptr(JFHCh,JfHH3E,JfHH3Ez,IUGvj,Z3XGkzi,VCoVV,VCoVVU,VCoVVUU),if (nargin<8),VCoVVUU=[0 0];end,g3E2mI=size(JfHH3E,1);g3E2mIN=size(JfHH3Ez,1);fQDmg=size(JfHH3Ez,2);if (~iscell(JFHCh)),JFHCh={JFHCh};end,[~,fQDmg1]=evalsfun(JFHCh{1},0,g3E2mI,g3E2mIN);fQDmg1k=sum(fQDmg1(:));if (nargin<6||isempty(VCoVV)),[tQRBbu,VCoVV]=mapdofbdr(JfHH3Ez,[],fQDmg1);else,fQDmg1k=size(VCoVV,1);tQRBbu=max(VCoVV(:));end,if (nargin<7||isempty(VCoVVU)),if (length(JFHCh)==1||strcmp(JFHCh{1},JFHCh{2})),iP5GR=fQDmg1k;tQRBb=tQRBbu;VCoVVU=VCoVV;else,[~,iP5GRj]=evalsfun(JFHCh{2},0,g3E2mI,g3E2mIN);iP5GR=sum(iP5GRj(:));[tQRBb,VCoVVU]=mapdofbdr(JfHH3Ez,[],iP5GRj);end,else,iP5GR=size(VCoVVU,1);tQRBb=max(VCoVVU(:));end,iP5GRj9=size(Z3XGkzi,2);cVoBQv=cell(1,iP5GRj9);cVoBQ=cell(1,iP5GRj9);CQcZU=ones(1,iP5GR);CQcZUh=ones(fQDmg1k,1)*[1:iP5GR];for CQcZUhH=1:iP5GRj9,fkvdGuk=Z3XGkzi{CQcZUhH};TeodBa5=VCoVV(:,fkvdGuk)'+VCoVVUU(2);TeodBa5=TeodBa5(:);TeodBa5=TeodBa5(:,CQcZU);cVoBQv{CQcZUhH}=TeodBa5(:);TeodBa5=VCoVVU(:,fkvdGuk)'+VCoVVUU(1);TeodBa5=TeodBa5(:,CQcZUh);cVoBQ{CQcZUhH}=TeodBa5(:);end,cVoBQv=cat(1,cVoBQv{:});cVoBQ=cat(1,cVoBQ{:});
function [ Z3XGkzi, qIzAl ]=asmgroups(qIzAlp,sind),Z3XGkzi={};qIzAlpb=max(sind);qIzAl=0;for fkvdGu=1:qIzAlpb,fkvdGuk=find(sind==fkvdGu);if (length(fkvdGuk)>qIzAlp),fkvdG=0;while (fkvdG<length(fkvdGuk)),VO7j7=fkvdG;fkvdG=min(fkvdG+qIzAlp,length(fkvdGuk));Z3XGkzi=[Z3XGkzi {fkvdGu;fkvdGuk(VO7j7+1:fkvdG)}];end,qIzAl=qIzAlp;else,Z3XGkzi=[Z3XGkzi {fkvdGu;fkvdGuk}];qIzAl=max(qIzAl,length(fkvdGuk));end,end,function [ T0VK7, T0VK7K, T0VK7Kb, NNX70l ]=evalsfun(NNX70l,g3E2m,g3E2mI,g3E2mIN,bMxdE,bMxdEl,bMxdElh,T0VK7),if (nargin==4),if (nargout==1),T0VK7=feval(NNX70l,g3E2m,g3E2mI,g3E2mIN);else,[T0VK7,T0VK7K,T0VK7Kb,NNX70l]=feval(NNX70l,g3E2m,g3E2mI,g3E2mIN);end,elseif (nargin>=7),if (nargin==7||isempty(T0VK7)),T0VK7(size(bMxdElh,1),1)=0;end,[T0VK7,T0VK7K,T0VK7Kb,NNX70l]=feval(NNX70l,g3E2m,g3E2mI,g3E2mIN,bMxdE,bMxdEl,bMxdElh,T0VK7);else,error('evalsfun: incorrect number of input arguments .'),end
