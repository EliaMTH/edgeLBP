%% Edge Local Binary Pattern (edgeLBP)
%
% a1=edgeLBPp(v,f,R,P,nrad,h) produces a matrix a1 (nv x R), each row
% represent the edgeLBP values of the vertex of the mesh. All the inputs
% are mandatory. The vertex and face list MUST be related to
% triangulation (with or without boundaries) which represent a SINGLE 
% manifold surface (with or without boundaries).
%
% [a1,data]=edgeLBPp(v,f,R,P,nrad,h) return also a matlab struct, used for
% debugging.
%
% If you use this are using this code, PLEASE cite us. You can find the
% bibtex in the same docs_and_readmes folder. 
%
% Copyright (c) Moscoso~Thompson Elia, 2018
%
% -------------------------------------------------------------------------
%
% INPUT and OUTPUT: 
%
% input:    - v    : vertex list, a matrix (nv x 3) or (3 x nv), where nv 
%                       is the number of vertices;
%           - f    : face list, a matrix (nf x 3) or (3 x nf), where nf is 
%                       the number of faces;
%           - R    : the size of the maximum radius;
%           - P    : the number of samples taken on each e-ring;
%           - nrad : the number of rings considered per vertex; 
%           - h    : the descriptor of the pattern, an arry with nv scalar
%                       entries, relative to the vertex list.
%
% output:   - a1   : the (nv x nrad) matrix, the j-row is the edgeLBP
%                       descriptor for the j-th vertex;
%           - data : the Matlab struct with some of the variables used. See
%                       output_data_readme.txt for details. Mostly used for
%                       debugging.
%

%%
function [a1,data]=edgeLBP(v,f,R,P,nrad,h)

% Vertex and face lists are trasposed in what we consider the standard
% shape. 
if size(v,2)~=3
    v=v';
end
if size(f,2)~=3
    f=f';
end
%

% Preliminary analysis on the mesh is now computed. All the results are
% stored in the prel structure. 
prel=eL_prel_tr(v,f,R,P,nrad);
prel.h=h;
tr=zeros(prel.nv,nrad);
a1=zeros(prel.nv,nrad);
vld=prel.vld;   % At this point of the code, all the vertices that are far 
                % enough from the boundary are valid. 
nv=prel.nv;
% 

% The edgeLBP values are now computed. 
% If a vertex is not admissible for any reason, the row that should
% report its edgeLBP values is replaced by a row of negative values, which
% give information on why that vertex is not admissible: 
% -1: the ring cannot be extracted; this is probably due to multiple
%       intestections between the sphere and the mesh;
% -2: the str array (see down in the code for the definition) cannot be
%       defined; this is probably due to the failing by the algorithm in
%       sorting the p_i. 
% -3: the vertex is too close to the boundary.
% The flag is 0 if the vertex is admissible.

parfor j=1:nv
    if vld(j)==1 
        [tr(j,:),a1(j,:),vld(j)]=f_loop(j,prel);
    else
        % FLAG3==Boundary vertex
        tr(j,:)=tr(j,:)-3;
        a1(j,:)=tr(j,:);
    end
end
%

% DEBUG output definition: all the variables are added to the data
% structure.
data.prel=prel;
data.prel.vld=vld;
data.tr=tr;
data.a1=a1;
%
end

%% ------------------------------------------------------------------------
function [ttr,ta1,vld]=f_loop(j,prel)
% Preallocating used variables
vld=1;
nrad=prel.nrad;
h=prel.h;
% This function computes the e-rings (in other words, computes the p_i)
[temp_ring,flag]=extract_rings(j,prel);  % Cell of rings
%

if flag==0
    % This sort the ring, without minding the starting point for now.
    temp_pi=cell(1,nrad);
    for r=1:nrad
        if flag~=-2
            [temp_pi{r},flag]=sorted_pi(temp_ring{r},prel);
        end
    end
end

if flag==0
    % FLAG = 0      ->    Admissible vertex
    
    % First, a starting point is give to each e-ring, then the s_i are
    % extracted from the p_i
    temp_pi=compute_des(temp_pi,prel); % compute h(p_i)
    str=get_str(temp_pi,prel,h(j));    % computing s_i
    
    % Initializing variables, then the transition and edgeLBP values are 
    % computed.
    ttr=zeros(1,prel.nrad); % This will contain the transition values
    ta1=zeros(1,prel.nrad); % This will contain the edgeLBP values.
    for r=1:nrad
        ttr(r)=get_tr(str{r});
        ta1(r)=sum(str{r});
    end
    %
    
    % if the flag is not 0, the respective flag is reported.
elseif flag==-1
    % FLAG = -1     ->     Problem in extracting rings
    vld=-1;
    ttr=zeros(1,prel.nrad)-1;
    ta1=ttr;
elseif flag==-2
    % FLAG = -2     ->     Problem in getting str
    vld=-2;
    ttr=zeros(1,prel.nrad)-2;
    ta1=ttr;
end


end
%% ------------------------------------------------------------------------
% Computing transitions, pretty straightforward.
function tr=get_tr(str)
nb=numel(str);
tr=0;
for j=2:nb
    if str(j-1)~=str(j)
        tr=tr+1;
    end
end
end

%% ------------------------------------------------------------------------
% Here the e-rings are computed, navigating the mesh vertices through the
% edges. In few words, a register of edges is created, which lists which
% edges have to be checked ('is this edge intersecting the sphere?') and 
% which have been already checked. The code is a bit convoluted, as I tried
% to make the minimum number of checks possibile. This is what takes the
% most computational time. 
% Also, if you are using this code to implement this function in other
% programming languages, this is probably not the optimal structure, from a
% performance point of view.

function [temp_ring,flag]=extract_rings(center,prel)

% The edge_register use the following labels:
% - 0: edges with no check required nor characterization
% - 1: edges to be checked
% - 2: edges already checked; any action on these vertices is dismissed.

% While I'm navigating the mesh, the label of the edges can change from 0 
% to 1 or from 1 to 2, based on the following code.

% Preallocating variables
edge_register=zeros(1,prel.ne);
edge_register(prel.ve{center})=1;   % The first edges to be checked are 
                                    % those that shares the vertex 
                                    % considered.
checking=find(edge_register==1,1);  % The index of the edge that is going 
                                    % to be checked.
dr=prel.R/prel.nrad;                % The size of the radii for the smaller
                                    % spheres.
flag=0;                             % Initializing the flag to 0.
temp_ring=cell(1,prel.nrad);
for j=1:prel.nrad                   % For each ring:
    temp_ring{j}.pi=[];             % p_i coodrinates are saved here
    temp_ring{j}.edg=[];            % which edges contain p_i is saved here
    temp_ring{j}.root=[];           % the value of the root t (see below)
end
ev=prel.e{1}(:,1:2);
pcenter=prel.v(center,:);           % Cartesian coordinates of the vertex 
                                    % considered.
todolist=[];                        % List of the edges that need to be 
                                    % checked.
%



while numel(checking)~=0
    % For this see "sphere_edge_inters.pdf". What follow is an
    % application of what is reported there.
    
    vtmp=ev(checking,1:2);
    v1=vtmp(1);
    v2=vtmp(2);
    pv1=prel.v(v1,:);
    pv2=prel.v(v2,:);
    
    for j=1:prel.nrad
        
        radius=dr*j;
        D=pv2-pcenter;
        c=pv1-pv2;
        alpha=norm(c)^2;
        beta2=sum(c.*D); %delta/4
        gamma=norm(D)^2-radius^2;
        d4=beta2^2-alpha*gamma;
        if d4<0
            t_root=[];
        else
            t_root=[(-beta2-sqrt(d4))/alpha,(-beta2+sqrt(d4))/alpha];
            ttemp=[];
            if t_root(1)>=0 && t_root(1)<=1
                ttemp=[ttemp,t_root(1)];
            end
            if t_root(2)>=0 && t_root(2)<=1
                ttemp=[ttemp,t_root(2)];
            end
            t_root=ttemp;
        end
        
        nroot=numel(t_root);
        
        switch nroot
            case 1
                temp_ring{j}.pi=[temp_ring{j}.pi;c*t_root+pv2];
                temp_ring{j}.edg=[temp_ring{j}.edg,checking];
                temp_ring{j}.root=[temp_ring{j}.root,t_root];
            case 2
                temp_ring{j}.pi=[temp_ring{j}.pi;c*t_root(1)+pv2;c*t_root(2)+pv2];
                temp_ring{j}.edg=[temp_ring{j}.edg,checking,checking];
                temp_ring{j}.root=[temp_ring{j}.root,t_root(1),t_root(2)];
        end
        for ll=1:nroot
            if t_root(ll)==0 || t_root(ll)==1
                flag=-1;
            end
        end
        
    end
    edge_register(checking)=2;
    if norm(pcenter-pv1)-prel.R<0 && norm(pcenter-pv2)-prel.R<0
        for jj=1:numel(prel.ee{checking})
            if edge_register(prel.ee{checking}(jj))==0
                edge_register(prel.ee{checking}(jj))=1;
            end
        end
    end
    
    if numel(todolist)==0
        todolist=find(edge_register==1);
    end
    if numel(todolist)>0
        checking=todolist(1);
        todolist(1)=[];
    else
        checking=[];
    end
    
    if flag==-1
        break;
    end
end

if exist('flag','var')==0
    flag=0;
end

end

%% ------------------------------------------------------------------------
% Getting sorted set of p_i. This is done two ways: in the standard
% situation the sorting is done using a face-to-face navigation of the
% mesh; if this is not possible, the closest point to the last sorted is
% considered as next. If after this procedure there are left out p_i, a new
% sorted ring is extracted from those. This is repeated untill no more
% p_i are left out. The bigger ring is then considered as the true e-ring.
%
% This process is the most efficient we found, in term of timing, 
% robustness, but it has a tricky requirement, checked before the sorting
% process starts.
% In more details, we check how many times i "visit" a  face, where a 
% "visit" is entering or leaving a the face. 
% In an optimal case, this number is 2. If a vertex of the triangulation is
% also a p_i, this number grow to 4. If both vertices of an edge are p_i,
% this number is 6, but this case is very rare and not considered. 

function [sorted_pile,flag]=sorted_pi(ring,prel)

ni=numel(ring.edg);
face_register=zeros(1,prel.nf);
flag=0;

for j=1:ni
    face_register(prel.ef(ring.edg(j),:))=face_register(prel.ef(ring.edg(j),:))+1;
end

recycle_flag=0;
% I'm sorting the p_i in order to have a vertex without root equal to 0 or
% 1.
while numel(find(ring.edg==ring.edg(1)))~=1 && recycle_flag<=ni
    ring.edg=ring.edg([2:end,1]);
    ring.pi=ring.pi([2:end,1],:);
    ring.root=ring.root([2:end,1]);
    recycle_flag=recycle_flag+1;
end

% This checks if it possible to sort the p_i.
if (max(face_register)==4 || max(face_register)==2) && numel(find(face_register==2))~=0
    
    raw_pile=prel.ef(ring.edg,:);
    cycle_flag=0;
    
    sorted_idx=zeros(1,ni);
    sorted_idx(1)=1;
    start=raw_pile(1,2);
    raw_pile(1,:)=[-1,-1];
    
    for j=2:ni
        [rr,cc]=find(raw_pile==start);
        switch numel(rr)
            case 1
                sorted_idx(j)=rr;
                start=raw_pile(rr,1+mod(cc,2));
                raw_pile(rr,:)=[-1,-1];
            case {2,3}
                nrr=numel(rr);
                pin=ring.pi(sorted_idx(j-1),:);
                d=zeros(1,nrr);
                for npos=1:nrr
                    d(npos)=norm(ring.pi(rr(npos),:)-pin);
                end
                closest=find(d==min(d));
                if numel(closest)>1
                    flag=-2;
                    break;
                end
                sorted_idx(j)=rr(closest);
                start=raw_pile(rr(closest),1+mod(cc(closest),2));
                raw_pile(rr(closest),:)=[-1,-1];
            otherwise
                Fcomp=sum(raw_pile(:,1)==-1);
                Scomp=sum(raw_pile(:,1)~=-1);
                sorted_idx=sorted_idx(1:j-1);
                if Fcomp<Scomp
                    cycle_flag=1;
                    break;
                end
                
                % ---------------------------------------------------------
                % Uncommenting this will show the sorted intersections.
                % First: all the intersections are shown, then they are
                % marked one by one as they are sorted.
                % ---------------------------------------------------------
                %                 figure
                %                 a=ring.pi;
                %                 hold on
                %                 for k=1:length(a)
                %                     plot3(a(k,1),a(k,2),a(k,3),'o','markersize',12);
                %                 end
                %                 a=ring.pi(sorted_idx,:);
                %                 hold on
                %                 for k=1:length(a)
                %                     plot3(a(k,1),a(k,2),a(k,3),'.','markersize',12);
                %                     pause(.3)
                %                 end
                %                 pause
                % ---------------------------------------------------------
                break;
        end
    end
    
    
    
    if cycle_flag==1
        
        sorted_idx=sorted_idx(1:j-1);
        temp_idx=1:ni;
        temp_idx(sorted_idx)=[];
        ring2.pi=ring.pi(temp_idx,:);
        ring2.edg=ring.edg(temp_idx);
        ring2.root=ring.root(temp_idx);
        [sorted_pile,flag]=sorted_pi(ring2,prel);
    end
    
    
    
else
    flag=-2;
    sorted_pile=-1;
end



if flag~=-2 && cycle_flag~=1
    sorted_pile.edg=ring.edg(sorted_idx);
    sorted_pile.pi=ring.pi(sorted_idx,:);
    sorted_pile.root=ring.root(sorted_idx);
end

end

%% ------------------------------------------------------------------------
% Getting str from sorted p_i
function str=get_str(temp_pi,prel,hcent)

pis=start_from_max(temp_pi); %Sorting the p_i
str=inter_sampl_confr(pis,prel,hcent); %Computing the s_i

end

%% ------------------------------------------------------------------------
% Compute descriptor h on the points p_i, pretty straightforward.
function pi=compute_des(pi,prel)

for jj=1:prel.nrad
    ni=numel(pi{jj}.root);
    pi{jj}.hpi=zeros(1,ni);
    for j=1:ni
        pi{jj}.hpi(j)=pi{jj}.root(j)*prel.h(prel.e{1}(pi{jj}.edg(j),1))+...
            (1-pi{jj}.root(j))*prel.h(prel.e{1}(pi{jj}.edg(j),2));
    end
    
end

end

%% ------------------------------------------------------------------------
% The starting point for the ring of the biggest radius is the one with the
% highest descriptor value, while for any other ring is the point closest
% to the initial point of the bigger ring.
% The implementation is straightforward.
function pis=start_from_max(pi)

nrad=numel(pi);

max_idx=find(pi{end}.hpi==max(pi{end}.hpi));
pis{nrad}.pi=[pi{end}.pi(max_idx:end,:);pi{end}.pi(1:max_idx-1,:)];
pis{nrad}.hpi=[pi{end}.hpi(max_idx:end),pi{end}.hpi(1:max_idx-1)];
ref_point=pi{end}.pi(max_idx,:);

for j=1:numel(pi)-1
    npi=size(pi{j}.pi,1);
    dists=zeros(1,npi);
    for jj=1:npi
        dists(jj)=norm(pi{j}.pi(jj,:)-ref_point);
    end
    
    max_idx=find(dists==min(dists));
    
    pis{j}.pi=[pi{j}.pi(max_idx:end,:);pi{j}.pi(1:max_idx-1,:)];
    pis{j}.hpi=[pi{j}.hpi(max_idx:end),pi{j}.hpi(1:max_idx-1)];
end

end

%% ------------------------------------------------------------------------
% This function output a P+1 array of either 0 or 1, as in the LBP
% algorithm. It takes a sorted e-ring with a starting point as input.
function str=inter_sampl_confr(pis,prel,hcent)

for j=1:prel.nrad
    samples{j}=inter_sampl(pis{j}.pi,pis{j}.hpi,prel.P);
end

for j=1:prel.nrad
    str{j}=zeros(1,prel.P);
    for k=1:prel.P
        if hcent<samples{j}(k)
            str{j}(k)=1;
        end
    end
end

end

%% ------------------------------------------------------------------------
% In this function, we extract P s_i samples from an e-ring. The linear
% interpolation iC between the p_i is considered, but it not explicitely
% computed. Instead, from a given starting point, iC is navigated until the
% lenght PP(iC)/P is reached (where PP is the perimeter): there a new 
% sample s_i is taken. This repeats until the whole Ic is navigated.
% Of course, some approximation are considered, in order to have a smoother
% and easier sampling.
% In more detail, we set s_1 as p_1. Then we "stack" multiple segment 
% (p_i,p_(i+1)) and we sum their length, until we reach PP(iC)/P. 
% If the sum is equal to that threshold (down to a 10^-3 * PP(iC)/P 
% precision), last p_i of the stack is a new sample. Then this cycle
% repeats, starting from the new sample. If by doing this I go over 
% PP(iC)/P, I create a sample in the segments that connects the two last
% p_i-s. 
% The following code is straightforward implementation of the method above. 

function samples=inter_sampl(pi,hpi,P)

PPiC=sum(sqrt(sum((pi-[pi(2:end,:);pi(1,:)]).^2,2)));

dr=PPiC/P;

np=numel(hpi);
samples=hpi(1);

k=1;%index on ring_des_rsl
j=1;
disc=0;
while k~=P
    jp=mod(j-1,np)+1;
    ja=mod(j,np)+1; 
    ddr=norm(pi(ja,:)-pi(jp,:));
    disc=disc+ddr; %<-distance covered
    if abs(disc-dr)<dr*(10^(-3))
        k=k+1;
        samples(k)=hpi(ja);
        j=ja;
        disc=0;
    elseif disc-dr>0
        k=k+1;
        lambd1=(disc-dr)/ddr;
        lambd2=1-lambd1;
        samples(k)=hpi(jp)*lambd1+hpi(ja)*lambd2;
        if jp==np
            hpi=[hpi(1:jp),samples(k)];
            pi=[pi(1:jp,:);pi(jp,:)*lambd1+pi(ja,:)*lambd2];
        else
            hpi=[hpi(1:jp),samples(k),hpi(ja:end)];
            pi=[pi(1:jp,:);pi(jp,:)*lambd1+pi(ja,:)*lambd2;pi(ja:end,:)];
        end
        np=np+1;
        j=jp+1;
        disc=0;
    else
        j=j+1;
    end
end


end

%% ------------------------------------------------------------------------
% Computing all the preliminary required relations between the entity of
% the triangulation, plus which elements are far enough from the boundary
% of the model.

function prel=eL_prel_tr(v,f,R,P,nrad)

% Computing face-face relation
fring=tri_adiac_cell(v,f); 

prel.P=P;
prel.R=R;
prel.v=v;
prel.f=f;
prel.nrad=nrad;
% Computing vertex-vertex relation
prel.vv=compute_vv_rel(f);

nf=length(f);
nv=length(v);

if size(f,1)==3
    f=f';
end
if size(v,1)==3
    v=v';
end

% Computing edges
ec=0;
for j=1:nf
    ec=ec+1;
    ed_l(ec)=norm(v(f(j,1),:)-v(f(j,2),:));
    edges(ec,:)=[sort([f(j,1),f(j,2)]),j,ed_l(ec)];
    
    ec=ec+1;
    ed_l(ec)=norm(v(f(j,2),:)-v(f(j,3),:));
    edges(ec,:)=[sort([f(j,2),f(j,3)]),j,ed_l(ec)];
    
    
    ec=ec+1;
    ed_l(ec)=norm(v(f(j,1),:)-v(f(j,3),:));
    edges(ec,:)=[sort([f(j,1),f(j,3)]),j,ed_l(ec)];
end

% Computing face-edge relation
[uniedg,wua,whereare]=unique(edges(:,1:2),'rows');
inx=zeros(1,nf);
prel.fe=zeros(nf,3);
for j=1:size(uniedg,1)
    prel.e{1}(j,1:3)=edges(wua(j),[1,2,4]);
    prel.e{2}{j}=edges(whereare==j,3);
    faceset=prel.e{2}{j};
    for jj=1:numel(faceset)
        inx(faceset(jj))=inx(faceset(jj))+1;
        prel.fe(faceset(jj),inx(faceset(jj)))=j;
    end
end

% Compute vertex-edge relation
prel.ve=compute_ve_rel(prel.e{1}(:,1:2));

% Computing boundary and internal vertices
valid=ones(nv,1);
nb=0;
ni=0;
board_vertex=[];
for j=1:nf
    if numel(fring{j}) < 3
        nb=nb+1;
        board_vertex((nb-1)*3+[1:3])=f(j,:);
        valid(f(j,:))=0;
    else
        ni=ni+1;
        inside_vertex((ni-1)*3+[1:3])=f(j,:);
    end
end
prel.bv=unique(board_vertex);
prel.iv=setdiff(unique(inside_vertex),prel.bv);
nb=numel(prel.bv);
ni=numel(prel.iv);


% Computing distance from board for all the inside vertex inside
vld=ones(1,nv);
if numel(prel.bv)~=0
    vld(prel.bv)=0;
    R2=R^2;
    for j=1:ni
        jj=1;
        while jj~=nb
            dist_temp=sum((v(prel.iv(j),:)-v(prel.bv(jj),:)).^2);
            if dist_temp < R2
                vld(prel.iv(j))=0;
                jj=nb;
            else
                jj=jj+1;
            end
        end
    end
end

prel.vld=vld;
prel.nv=nv;
prel.ni=ni;
prel.nb=nb;
prel.nf=nf;

% Computing other triangulation relations
prel.ve=compute_ve_rel(prel.e{1}(:,1:2));
prel.ne=numel(prel.e{2});
prel.ef=compute_ef_rel(prel.fe);
prel.ee=compute_ee_rel(prel.e{1}(:,1:2),prel.ve);
end

%% ------------------------------------------------------------------------
% Functions for triangulation relations
function ee=compute_ee_rel(edge,ve)
for j=1:length(edge)
    ee{j}=horzcat(ve{edge(j,1)},ve{edge(j,2)});
end
end
%
function ef=compute_ef_rel(fe)
ef=zeros(max(max(fe)),2);
for j=1:max(max(fe))
    [r,~]=find(fe==j);
    if numel(r)==1
        ef(j,:)=[r,-1];
    else
        ef(j,:)=r';
    end
end
end
%
function [ad_tri]=tri_adiac_cell(mat_ver,mat_tri)
num_ver=size(mat_ver,1); num_tri=size(mat_tri,1);

succ=[2 3 1]; t=[1:num_tri];

t_t=spalloc(num_ver,num_ver,3*num_tri);
for i=1:3
    t_t=t_t+sparse(mat_tri(:,i),mat_tri(:,succ(i)),t,num_ver,num_ver,num_tri);
end

for i=1:3
    index_edge=sub2ind(size(t_t),mat_tri(:,succ(i)),mat_tri(:,i));
    mat_tri(:,i+3)=t_t(index_edge);
end

mat_tri=mat_tri(:,4:6);
ad_tri=cell(1,size(mat_tri,1));

for i=1:size(mat_tri,1)
    bb=[];
    cc=0;
    for ii=1:3
        if mat_tri(i,ii)~=0
            cc=cc+1;
            bb(1,cc)=mat_tri(i,ii);
        end
    end
    cc=0;
    ad_tri{1,i}=bb;
end

end
%
function ve=compute_ve_rel(edge)

nv=numel(min(edge(:)):max(edge(:)));
ne=length(edge);
ve=cell(nv,1);

for j=1:ne
    ve{edge(j,1)}=[ve{edge(j,1)},j];   
    ve{edge(j,2)}=[ve{edge(j,2)},j];  
end

for j=1:nv
    ve{j}=unique(ve{j});
end
end
%
function vv_rel=compute_vv_rel(f)

nf=size(f,1);
nv=max(max(f));
index1=zeros(nv,1);

for j=1:nf
    for i=1:3
        index1(f(j,i))=index1(f(j,i))+1;
        vv_rel{f(j,i)}(index1(f(j,i)))=f(j,mod(i,3)+1);
        
        index1(f(j,i))=index1(f(j,i))+1;
        vv_rel{f(j,i)}(index1(f(j,i)))=f(j,mod(i+1,3)+1);
    end
end

for j=1:numel(vv_rel)
    vv_rel{j}=unique(vv_rel{j});
end


end
