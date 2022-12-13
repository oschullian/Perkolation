clear all
close all


rng(2394)

np=100
pall=linspace(0.4,0.6,np);
p2=linspace(0,1,np)
ns0=5
nsall=[ns0,ns0,ns0]

for i=1:length(nsall)
    ns=nsall(i)
    data(i).Q=zeros(ns^2*2+1,ns^2*2)
    data(i).percol_occur=zeros(1,ns^2*2+1)
    data(i).percolatedclsize=zeros(1,ns^2*2+1)
    data(i).percolatedclsizerun=zeros(1,ns^2*2+1)
    data(i).run=0
    data(i).nofiniteclusters=zeros(1,ns^2*2+1)
end

Pperc=zeros(length(nsall),np);
Psize=zeros(length(nsall),np);
Qdist=zeros(length(nsall),max(nsall)^2*2);
K=zeros(length(nsall),np);


for irun=1:1000
    irun
    for is=1:length(nsall)
        data(is).run=data(is).run+1;
        ns=nsall(is);
        %percol_occur=zeros(1,ns^2*2+1);
        allbonds=generatebonds(ns);
        occ=zeros(2,ns,ns);
        clname=nan(2,ns,ns);
        il=0;
        iu=0;
        Qn=zeros(ns^2*2+1,ns^2*2);
        
        for in=1:ns^2*2
            
            il=iu+1;
            iu=in;
            
            [clname,occ]=findclusters(occ,clname,allbonds,il,iu);
            
            [clname_periodic,offsetx,offsety,pernames]=findclusters_periodic(clname,occ);
            
            plotclusters_percolation(clname_periodic,occ,1,pernames)
            
            if ~isempty(pernames)
                data(is).percol_occur(in+1)=data(is).percol_occur(in+1)+length(pernames);
                for i=1:length(pernames)
                    data(is).percolatedclsize(in+1)=data(is).percolatedclsize(in+1)+sum(clname_periodic(:)==pernames(i));
                end
                data(is).percolatedclsizerun(in+1)=data(is).percolatedclsizerun(in+1)+1;
            end
            
            clnames_un=unique(clname_periodic(~isnan(clname_periodic)));
            if ~isempty(pernames)
                for i=1:length(pernames)
                    clnames_un=clnames_un(clnames_un~=pernames(i));
                end
            end
            for i=1:length(clnames_un)
                ind=sum(clname_periodic(:)==clnames_un(i));
                Qn(in+1,ind)=Qn(in+1,ind)+1;
            end
            
            data(is).nofiniteclusters(in+1)=data(is).nofiniteclusters(in+1)+length(clnames_un);
            
        end
        
        
        
        Qnsum=sum(Qn,2);
        for i=1:ns^2*2+1
            if(Qnsum(i)~=0)
                Qn(i,:)=Qn(i,:)/Qnsum(i);
            end
        end
        data(is).Q=data(is).Q+Qn;
 
    end
    save(sprintf('data%i.mat',ns0),'data')
end

% Die Funktion weisst eine Kante einem Cluster zu
function [clname,occ]=findclusters(occ,clname,allbonds,il,iu)


ns=size(occ,2);

% running index
clnp=max(clname(:));
if isnan(clnp)
    clnp=0;
end
% durchgehen aller Kanten um zu schauen ob sie mit den benachbarten Kanten
% verbunden sind


for ib=il:iu
    
    id=allbonds(ib,1);
    ix=allbonds(ib,2);
    iy=allbonds(ib,3);
    occ(id,ix,iy)=1;
    
    % zunächst betrachten wir die horizontalen Kanten. Jede dieser
    % Kanten hat drei benachbarte Kanten an den beiden Knoten. Hier
    % wird getestet ob diese offen sind.
    if id==1
        % teste of Kante offen ist
        if occ(1,ix,iy)==0
            error('bond is closed')
        end
        % falls die Kante noch keinem Cluster zugeordnet ist, ordne sie
        % einem neuen Cluster zu
        if isnan(clname(1,ix,iy))
            clnp=clnp+1;
            clname(1,ix,iy)=clnp;
        end
        
        % nun werden alle benachbarten Kanten betrachtet und es wird
        % getested ob diese offen sind
        
        % linker Knoten, Kante nach links
        if ix>=2
            % falls die benachbarte Kante offen ist, so gib ihr und
            % allen Kanten im selben Cluster den Knotennamen der gerade
            % betrachten Kante
            if occ(1,ix-1,iy)==1
                if isnan(clname(1,ix-1,iy))
                    clname(1,ix-1,iy)=clname(1,ix,iy);
                else
                    clname(clname(:)==clname(1,ix-1,iy))=clname(1,ix,iy);
                end
            end
        end
        
        % linker Knoten, Kante nach oben
        if occ(2,ix,iy)==1
            if isnan(clname(2,ix,iy))
                clname(2,ix,iy)=clname(1,ix,iy);
            else
                clname(clname(:)==clname(2,ix,iy))=clname(1,ix,iy);
            end
        end
        
        % linker Knoten, Kante nach unten
        if iy>=2
            if occ(2,ix,iy-1)==1
                if isnan(clname(2,ix,iy-1))
                    clname(2,ix,iy-1)=clname(1,ix,iy);
                else
                    clname(clname(:)==clname(2,ix,iy-1))=clname(1,ix,iy);
                end
            end
        end
        
        % rechter Knoten, Kante nach unten
        if ix<ns
            if iy>=2
                if occ(2,ix+1,iy-1)==1
                    if isnan(clname(2,ix+1,iy-1))
                        clname(2,ix+1,iy-1)=clname(1,ix,iy);
                    else
                        clname(clname(:)==clname(2,ix+1,iy-1))=clname(1,ix,iy);
                    end
                end
            end
        end
        
        % rechter Knoten, Kante nach rechts
        if ix<ns
            if occ(1,ix+1,iy)==1
                if isnan(clname(1,ix+1,iy))
                    clname(1,ix+1,iy)=clname(1,ix,iy);
                else
                    clname(clname(:)==clname(1,ix+1,iy))=clname(1,ix,iy);
                end
            end
        end
        
        % rechter Knoten, Kante nach oben
        if ix<ns
            if occ(2,ix+1,iy)==1
                if isnan(clname(2,ix+1,iy))
                    clname(2,ix+1,iy)=clname(1,ix,iy);
                else
                    clname(clname(:)==clname(2,ix+1,iy))=clname(1,ix,iy);
                end
            end
        end
        
        
    else
        
        % teste of Kante offen ist
        if occ(2,ix,iy)==0
            error('bond is closed')
        end
        
        
        if isnan(clname(2,ix,iy))
            clnp=clnp+1;
            clname(2,ix,iy)=clnp;
        end
        
        
        % unterer Knoten, Kante nach rechts
        if occ(1,ix,iy)==1
            if isnan(clname(1,ix,iy))
                clname(1,ix,iy)=clname(2,ix,iy);
            else
                clname(clname(:)==clname(1,ix,iy))=clname(2,ix,iy);
            end
        end
        
        % unterer Knoten, Kante nach links
        if ix>=2
            if occ(1,ix-1,iy)==1
                if isnan(clname(1,ix-1,iy))
                    clname(1,ix-1,iy)=clname(2,ix,iy);
                else
                    clname(clname(:)==clname(1,ix-1,iy))=clname(2,ix,iy);
                end
            end
        end
        
        % unterer Knoten, Kante nach unten
        if iy>=2
            if occ(2,ix,iy-1)==1
                if isnan(clname(2,ix,iy-1))
                    clname(2,ix,iy-1)=clname(2,ix,iy);
                else
                    clname(clname(:)==clname(2,ix,iy-1))=clname(2,ix,iy);
                end
            end
        end
        
        
        % oberer Knoten, Kante nach links
        if iy<ns
            if ix>=2
                if occ(1,ix-1,iy+1)==1
                    if isnan(clname(1,ix-1,iy+1))
                        clname(1,ix-1,iy+1)=clname(2,ix,iy);
                    else
                        clname(clname(:)==clname(1,ix-1,iy+1))=clname(2,ix,iy);
                    end
                end
            end
        end
        
        % oberer Knoten, Kante nach oben
        if iy<ns
            if occ(2,ix,iy+1)==1
                if isnan(clname(2,ix,iy+1))
                    clname(2,ix,iy+1)=clname(2,ix,iy);
                else
                    clname(clname(:)==clname(2,ix,iy+1))=clname(2,ix,iy);
                end
            end
        end
        
        % oberer Knoten, Kante nach rechts
        if iy<ns
            if occ(1,ix,iy+1)==1
                if isnan(clname(1,ix,iy+1))
                    clname(1,ix,iy+1)=clname(2,ix,iy);
                else
                    clname(clname(:)==clname(1,ix,iy+1))=clname(2,ix,iy);
                end
            end
        end
        
        
    end
    
end


clnames_un=unique(clname(~isnan(clname)));

for i=1:length(clnames_un)
    clname(clname(:)==clnames_un(i))=i;
end



end




% this function generates a table with the direction (1 for horizontal, 2
% for vertical) and the poisition of the lattice site to the left or to the
% right (depending the direction of the bond). The last columns is a random
% uniformly distributed number (to decide if a bond is created or not)
function allbonds=generatebonds(ns)
ixall=zeros(2,ns,ns);
iyall=zeros(2,ns,ns);
dirall=zeros(2,ns,ns);
dirall(1,:,:)=1;
dirall(2,:,:)=2;
for ix=1:ns
    ixall(:,ix,:)=ix;
    iyall(:,:,ix)=ix;
end
allbonds=[dirall(:),ixall(:),iyall(:)];
[a,b]=sort(rand(ns^2*2,1));
allbonds=allbonds(b,:);
end


% finde cluster, die über die periodischen Randbedingungen gehen
% Falls eine Kante zwei vorher unabhängige Cluster über die Boxgrenze hin
% verknüpft, wird eines der Cluster so verschoben (mit einem
% Verschiebungsvektor), sodass das Cluster
% vervollständigt wird. Falls eine Bindung zwei Cluster von derselben
% Bezeichnung verbindet wird getestet, ob der Verschiebungsvektor genau
% einer Verschiebung um eine Box in die entsprechende Richtung entspricht.
% Falls ja, ist das Cluster endlich (es geht nur über die Box), falls
% nicht, ist es ein unendliches Cluster, d.h. Perkolation findet statt
function [clname_periodic,offsetx,offsety,pernames]=findclusters_periodic(clname,occ)
ns=size(clname,2);
offsetx=zeros(2,ns,ns);
offsety=zeros(2,ns,ns);
clname_periodic=clname;

% name der perkolierenden Cluster
pernames=[];

% betrachte die Kanten, die über die obere periodische Boxgrenze gehen
for ix=1:ns
    
    % betrachte nur Kanten welche offen sind
    if occ(2,ix,ns)==1
        
        % falls die Kante noch nicht einer Box zugeordnet sind, weise sie
        % der zentralen Box zu
        if isnan(offsetx(2,ix,ns))
            offsetx(clname(:)==clname(2,ix,ns))=0;
            offsety(clname(:)==clname(2,ix,ns))=0;
        end
        
        
        % oberer Knoten, Kante nach oben
        if occ(2,ix,1)==1
            % falls es sich um zwei unterschiedliche Cluster handelt,
            % können wir das gesamt Cluster verschieben ohne die eben
            % betrachtete Kante zu verschieben
            if clname_periodic(2,ix,1)~=clname_periodic(2,ix,ns)
                
                offsetx(clname_periodic(:)==clname_periodic(2,ix,1))=...
                    offsetx(clname_periodic(:)==clname_periodic(2,ix,1))+...
                    offsetx(2,ix,ns)-offsetx(2,ix,1);
                offsety(clname_periodic(:)==clname_periodic(2,ix,1))=offsety(clname_periodic(:)==clname_periodic(2,ix,1))+...
                    offsety(2,ix,ns)+1-offsety(2,ix,1);
                
                % verknüpfe die beiden Cluster, indem alle Kanten gleich
                % bezeichnet werden
                pernames(pernames==clname_periodic(2,ix,1))=clname_periodic(2,ix,ns);
                clname_periodic(clname_periodic(:)==clname_periodic(2,ix,1))=...
                    clname_periodic(2,ix,ns);
            else
                % es handelt sich bereits um dasselbe Cluster, entsprechen
                % die Verschiebungsvektoren bereits ver korrekten
                % Verschiebung oder handelt es sich um Perkolation?
                if not(and((offsety(2,ix,1)-offsety(2,ix,ns))==1,(offsetx(2,ix,1)-offsetx(2,ix,ns))==0))
                    % speichere den Namen des Clusters, welches perkoliert
                    pernames=unique([pernames,clname_periodic(2,ix,ns)]);
                end
            end
        end
        
        
        % oberer Knoten, Kante nach rechts
        if occ(1,ix,1)==1
            if clname_periodic(1,ix,1)~=clname_periodic(2,ix,ns)
                offsetx(clname_periodic(:)==clname_periodic(1,ix,1))=...
                    offsetx(clname_periodic(:)==clname_periodic(1,ix,1))+...
                    offsetx(2,ix,ns)-offsetx(1,ix,1);
                
                offsety(clname_periodic(:)==clname_periodic(1,ix,1))=...
                    offsety(clname_periodic(:)==clname_periodic(1,ix,1))+...
                    offsety(2,ix,ns)+1-offsety(1,ix,1);
                
                pernames(pernames==clname_periodic(1,ix,1))=clname_periodic(2,ix,ns);
                clname_periodic(clname_periodic(:)==clname_periodic(1,ix,1))=...
                    clname_periodic(2,ix,ns);
            else
                if not(and((offsety(1,ix,1)-offsety(2,ix,ns))==1,(offsetx(1,ix,1)-offsetx(2,ix,ns))==0))
                    pernames=unique([pernames,clname_periodic(2,ix,ns)]);
                end
            end
        end
        
        % oberer Knoten, Kante nach links
        
        % betrachte, die Fallunterscheidung für den Fall, dass die offene
        % Verbindung über zwei Boxen geht, d.h. in der linken oberen Ecke
        if ix>1
            if occ(1,ix-1,1)==1
                if clname_periodic(1,ix-1,1)~=clname_periodic(2,ix,ns)
                    offsetx(clname_periodic(:)==clname_periodic(1,ix-1,1))=...
                        offsetx(clname_periodic(:)==clname_periodic(1,ix-1,1))+...
                        offsetx(2,ix,ns)-offsetx(1,ix-1,1);
                    
                    offsety(clname_periodic(:)==clname_periodic(1,ix-1,1))=...
                        offsety(clname_periodic(:)==clname_periodic(1,ix-1,1))+...
                        offsety(2,ix,ns)+1-offsety(1,ix-1,1);
                    
                    pernames(pernames==clname_periodic(1,ix-1,1))=clname_periodic(2,ix,ns);
                    clname_periodic(clname_periodic(:)==clname_periodic(1,ix-1,1))=...
                        clname_periodic(2,ix,ns);
                else
                    if not(and((offsety(1,ix-1,1)-offsety(2,ix,ns))==1,(offsetx(1,ix-1,1)-offsetx(2,ix,ns))==0))
                        pernames=unique([pernames,clname_periodic(2,ix,ns)]);
                    end
                end
            end
        else
            if occ(1,ns,1)==1
                if clname_periodic(1,ns,1)~=clname_periodic(2,ix,ns)
                    offsetx(clname_periodic(:)==clname_periodic(1,ns,1))=...
                        offsetx(clname_periodic(:)==clname_periodic(1,ns,1))+...
                        offsetx(2,ix,ns)-offsetx(1,ns,1)-1;
                    
                    offsety(clname_periodic(:)==clname_periodic(1,ns,1))=...
                        offsety(clname_periodic(:)==clname_periodic(1,ns,1))+...
                        offsety(2,ix,ns)+1-offsety(1,ns,1);
                    
                    pernames(pernames==clname_periodic(1,ns,1))=clname_periodic(2,ix,ns);
                    clname_periodic(clname_periodic(:)==clname_periodic(1,ns,1))=...
                        clname_periodic(2,ix,ns);
                else
                    if not(and((offsety(1,ns,1)-offsety(2,ix,ns))==1,(offsetx(1,ns,1)-offsetx(2,ix,ns))==-1))
                        pernames=unique([pernames,clname_periodic(2,ix,ns)]);
                    end
                end
            end
        end
    end
end


% betrachte die Kanten, die über die rechte periodische Boxgrenze gehen
for iy=1:ns
    if occ(1,ns,iy)==1
        
        if isnan(offsetx(1,ns,iy))
            offsetx(clname(:)==clname(1,ns,iy))=0;
            offsety(clname(:)==clname(1,ns,iy))=0;
        end
        
        % rechter Knoten, Kante nach oben
        if occ(2,1,iy)==1
            if clname_periodic(2,1,iy)~=clname_periodic(1,ns,iy)
                offsetx(clname_periodic(:)==clname_periodic(2,1,iy))=...
                    offsetx(clname_periodic(:)==clname_periodic(2,1,iy))+...
                    offsetx(1,ns,iy)-offsetx(2,1,iy)+1;
                offsety(clname_periodic(:)==clname_periodic(2,1,iy))=offsety(clname_periodic(:)==clname_periodic(2,1,iy))+...
                    offsety(1,ns,iy)-offsety(2,1,iy);
                
                pernames(pernames==clname_periodic(2,1,iy))=clname_periodic(1,ns,iy);
                clname_periodic(clname_periodic(:)==clname_periodic(2,1,iy))=...
                    clname_periodic(1,ns,iy);
            else
                if not(and((offsety(2,1,iy)-offsety(1,ns,iy))==0,(offsetx(2,1,iy)-offsetx(1,ns,iy))==1))
                    pernames=unique([pernames,clname_periodic(1,ns,iy)]);
                end
            end
        end
        
        
        
        % rechter Knoten, Kante nach rechts
        if occ(1,1,iy)==1
            if clname_periodic(1,1,iy)~=clname_periodic(1,ns,iy)
                
                offsetx(clname_periodic(:)==clname_periodic(1,1,iy))=...
                    offsetx(clname_periodic(:)==clname_periodic(1,1,iy))+...
                    offsetx(1,ns,iy)-offsetx(1,1,iy)+1;
                
                offsety(clname_periodic(:)==clname_periodic(1,1,iy))=...
                    offsety(clname_periodic(:)==clname_periodic(1,1,iy))+...
                    offsety(1,ns,iy)-offsety(1,1,iy);
                
                pernames(pernames==clname_periodic(1,1,iy))=clname_periodic(1,ns,iy);
                clname_periodic(clname_periodic(:)==clname_periodic(1,1,iy))=...
                    clname_periodic(1,ns,iy);
            else
                if and((offsety(1,1,iy)-offsety(1,ns,iy))==0,(offsetx(1,1,iy)-offsetx(1,ns,iy))==1)
                else
                    pernames=unique([pernames,clname_periodic(1,ns,iy)]);
                end
            end
        end
        
        % rechter Knoten, Kante nach unten
        if iy>=2
            if occ(2,1,iy-1)==1
                if clname_periodic(2,1,iy-1)~=clname_periodic(1,ns,iy)
                    
                    offsetx(clname_periodic(:)==clname_periodic(2,1,iy-1))=...
                        offsetx(clname_periodic(:)==clname_periodic(2,1,iy-1))+...
                        offsetx(1,ns,iy)-offsetx(2,1,iy-1)+1;
                    offsety(clname_periodic(:)==clname_periodic(2,1,iy-1))=...
                        offsety(clname_periodic(:)==clname_periodic(2,1,iy-1))+...
                        offsety(1,ns,iy)-offsety(2,1,iy-1);
                    
                    pernames(pernames==clname_periodic(2,1,iy-1))=clname_periodic(1,ns,iy);
                    clname_periodic(clname_periodic(:)==clname_periodic(2,1,iy-1))=...
                        clname_periodic(1,ns,iy);
                else
                    if not(and((offsety(2,1,iy-1)-offsety(1,ns,iy))==0,(offsetx(2,1,iy-1)-offsetx(1,ns,iy))==1))
                        pernames=unique([pernames,clname_periodic(1,ns,iy)]);
                    end
                end
            end
        else
            if occ(2,1,ns)==1
                if clname_periodic(2,1,ns)~=clname_periodic(1,ns,iy)
                    
                    offsetx(clname_periodic(:)==clname_periodic(2,1,ns))=...
                        offsetx(clname_periodic(:)==clname_periodic(2,1,iy-1))+...
                        offsetx(1,ns,iy)-offsetx(2,1,ns)+1;
                    offsety(clname_periodic(:)==clname_periodic(2,1,ns))=...
                        offsety(clname_periodic(:)==clname_periodic(2,1,iy-1))+...
                        offsety(1,ns,iy)-offsety(2,1,ns)-1;
                    
                    pernames(pernames==clname_periodic(2,1,ns))=clname_periodic(1,ns,iy);
                    clname_periodic(clname_periodic(:)==clname_periodic(2,1,ns))=...
                        clname_periodic(1,ns,iy);
                else
                    if not(and((offsety(2,1,ns)-offsety(1,ns,iy))==-1,(offsetx(2,1,ns)-offsetx(1,ns,iy))==1))
                        pernames=unique([pernames,clname_periodic(1,ns,iy)]);
                    end
                end
            end
        end
    end
end

pernames=unique(pernames);
clnames_un=unique(clname_periodic(~isnan(clname_periodic)));
for i=1:length(clnames_un)
    pernames(pernames==clnames_un(i))=i;
    clname_periodic(clname_periodic(:)==clnames_un(i))=i;
end


end




% this function plots the grid
function plotclusters(clname,occ,fign)
% plot just occupied bonds
t=linspace(0,1,6);
colm=[0.1,1,1
    0,0.8,0
    1,1,0
    1,0,0
    1,0,1
    0,0,0.8];

tint=((0:max(clname(:)))/(max(clname(:))))';
col=[interp1(t,colm(:,1),tint),...
    interp1(t,colm(:,2),tint),...
    interp1(t,colm(:,3),tint)];
figure(fign)
hold off
for ix=1:size(clname,2)
    for iy=1:size(clname,2)
        if occ(1,ix,iy)==1
            plot([ix-1,ix],[iy-1,iy-1],'Linewidth',2,'Color',col(clname(1,ix,iy),:))
            hold on
        end
        if occ(2,ix,iy)==1
            plot([ix-1,ix-1],[iy-1,iy],'Linewidth',2,'Color',col(clname(2,ix,iy),:))
            hold on
        end
    end
end
drawnow
end


function plotclusters_percolation(clname,occ,fign,pernames)
% plot just occupied bonds


t=linspace(0,1,6);
colm=[0.1,1,1
    0,0.8,0
    1,1,0
    1,0,0
    1,0,1
    0,0,0.8];

tint=((0:max(clname(:)))/(max(clname(:))))';
col=[interp1(t,colm(:,1),tint),...
    interp1(t,colm(:,2),tint),...
    interp1(t,colm(:,3),tint)];

figure(fign)
hold off

for ix=1:size(clname,2)
    for iy=1:size(clname,2)
        
        if occ(1,ix,iy)==1
            if any(clname(1,ix,iy)==pernames)
                lw=3;
            else
                lw=2;
            end
            plot([ix-1,ix],[iy-1,iy-1],'Linewidth',lw,'Color',col(clname(1,ix,iy),:))
            hold on
        end
        
        if occ(2,ix,iy)==1
            if any(clname(2,ix,iy)==pernames)
                lw=3;
            else
                lw=2;
            end
            plot([ix-1,ix-1],[iy-1,iy],'Linewidth',lw,'Color',col(clname(2,ix,iy),:))
            hold on
        end
    end
end
drawnow
end