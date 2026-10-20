%% EXPERIMENTO 4: Givens hibrido x FastICA x JADE x SOBI
% MATLAB R2016b+. Requer implementacoes externas: fastica.m, jadeR.m, sobi.m
% Instale as implementacoes de referencia no MATLAB path ANTES de executar.
% Nao substitui algoritmos ausentes por aproximacoes nao equivalentes.
clear; clc; rng(20261020,'twister');
N=6000; numTrials=20; numStarts=20; delays=[1 2 3 5 8]; rho=10;
methods={'GivensHibrido','FastICA','JADE','SOBI'};
scenarios={'A_distintos','B_temporal_igual','C_quase_gauss','D_complementar'};
required={'fastica','jadeR','sobi'};
for i=1:numel(required)
    if exist(required{i},'file')~=2
        error('Funcao %s.m ausente no path. Instale a implementacao de referencia e execute novamente.',required{i});
    end
end
% Conferir convencoes das implementacoes instaladas (ver adaptadores abaixo).
opts=optimset('Display','off','MaxFunEvals',2500,'MaxIter',1500,'TolX',1e-8,'TolFun',1e-10);
% results: sc,trial,method,[MeanCorr,Amari,tempoSegundos]
results=nan(4,numTrials,4,3);
for sc=1:4
 fprintf('\n=== %s ===\n',scenarios{sc});
 for trial=1:numTrials
    s=makeSources(sc,N);
    [U,~]=qr(randn(3)); [V,~]=qr(randn(3));
    A=U*diag([1 1.4 2])*V';
    if det(A)<0, A(:,1)=-A(:,1); end
    x=A*s; x=x-mean(x,2);
    [E,L]=eig((x*x')/N);
    whitening=diag(1./sqrt(max(diag(L),eps)))*E';
    z=whitening*x;
    C=zeros(3,3,numel(delays));
    for k=1:numel(delays)
      lag=delays(k); R=z(:,lag+1:end)*z(:,1:end-lag)'/(N-lag);
      C(:,:,k)=(R+R')/2;
    end
    initialAngles=2*pi*rand(numStarts,3);
    for m=1:4
      tic;
      try
        switch m
          case 1
            best=-inf; B=eye(3);
            for st=1:numStarts
              fun=@(th)-hybridObjective(givensQ(th),z,C,rho);
              th=fminsearch(fun,initialAngles(st,:),opts);
              val=-fun(th);
              if isfinite(val) && val>best
                best=val; B=givensQ(th);
              end
            end
          case 2
            % fastica(icasig input) pode rebranquear internamente.
            % 'whiteSig' pede usar os dados z ja branqueados em versoes padrao.
            % Verifique a assinatura da versao instalada.
            [~,~,B]=fastica(z,'approach','symm','g','tanh',...
                'numOfIC',3,'verbose','off','displayMode','off');
          case 3
            % jadeR retorna matriz separadora (n x n), aplicada a z.
            B=jadeR(z,3);
          case 4
            % Implementacao instalada: [H,S,D]=sobi(X,n,p).
            % H e matriz de mistura estimada; usar inversa para separacao.
            % p e numero de defasagens (escalar), nao vetor de atrasos.
            [H,~,~]=sobi(z,3,max(delays));
            B=pinv(H);
        end
        elapsed=toc;
        if ~isequal(size(B),[3 3]) || any(~isfinite(B(:)))
           error('Matriz de separacao invalida: verifique assinatura da biblioteca.');
        end
        y=B*z;
        [mc,am]=metrics(y,s,B*whitening*A);
        results(sc,trial,m,:)=[mc am elapsed];
      catch ME
        elapsed=toc;
        warning('%s, trial %d, %s: %s',scenarios{sc},trial,methods{m},ME.message);
        results(sc,trial,m,3)=elapsed;
      end
    end
    fprintf('  %02d/%02d\n',trial,numTrials);
 end
end
fid=fopen('ica_givens_experimento4.csv','w');
fprintf(fid,'cenario,trial,metodo,MeanCorr,Amari,tempo_s\n');
for sc=1:4
 for t=1:numTrials
  for m=1:4
   v=reshape(results(sc,t,m,:),1,3);
   fprintf(fid,'%s,%d,%s,%.12g,%.12g,%.12g\n',scenarios{sc},t,methods{m},v);
  end
 end
end
fclose(fid);
fid=fopen('ica_givens_experimento4_resumo.csv','w');
fprintf(fid,'cenario,metodo,n_validos,MeanCorr_media,MeanCorr_dp,Amari_media,Amari_dp,tempo_media_s,tempo_dp_s\n');
for sc=1:4
 fprintf('\n%s\n',scenarios{sc});
 for m=1:4
  a=reshape(results(sc,:,m,1),[],1); b=reshape(results(sc,:,m,2),[],1);
  c=reshape(results(sc,:,m,3),[],1); valid=isfinite(a)&isfinite(b)&isfinite(c);
  if any(valid)
   av=mean(a(valid)); ad=std(a(valid)); bv=mean(b(valid)); bd=std(b(valid));
   cv=mean(c(valid)); cd=std(c(valid));
  else
   av=NaN;ad=NaN;bv=NaN;bd=NaN;cv=NaN;cd=NaN;
  end
  fprintf('%-15s n=%d  Corr=%.5f  Amari=%.5f  Tempo=%.4fs\n',methods{m},sum(valid),av,bv,cv);
  fprintf(fid,'%s,%s,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',...
   scenarios{sc},methods{m},sum(valid),av,ad,bv,bd,cv,cd);
 end
end
fclose(fid);
save('ica_givens_experimento4.mat','results','methods','scenarios','N','numTrials','numStarts','delays','rho');
fprintf('\nArquivos CSV e MAT salvos. Conferir avisos de adaptadores.\n');

function s=makeSources(sc,N)
 lap=(-log(rand(1,N))+log(rand(1,N)))/sqrt(2);
 uni=sqrt(3)*(2*rand(1,N)-1);
 rad=2*(rand(1,N)>.5)-1;
 switch sc
 case 1
  u=[lap;uni;rad];
  s=[filter(1,[1 -.65],u(1,:));filter(1,[1 .55],u(2,:));filter([1 .7],1,u(3,:))];
 case 2
  u=[lap;uni;rad];s=zeros(3,N);
  for i=1:3, s(i,:)=filter(1,[1 -.65],u(i,:));end
 case 3
  r=.05;u=sqrt(1-r^2)*randn(3,N)+r*[lap;uni;rad];
  s=[filter(1,[1 -.65],u(1,:));filter(1,[1 .55],u(2,:));filter([1 .7],1,u(3,:))];
 case 4
  u=[lap;uni;rad];
  s=[filter(1,[1 -.65],u(1,:));filter(1,[1 -.65],u(2,:));filter(1,[1 .55],u(3,:))];
 end
 s=s-mean(s,2);s=bsxfun(@rdivide,s,std(s,0,2));
end
function J=hybridObjective(Q,z,C,rho)
 y=Q*z;y=y-mean(y,2);v=mean(y.^2,2);
 kurt=mean(y.^4,2)./(v.^2)-3;
 K=sum(kurt.^2);D=0;
 for k=1:size(C,3)
  B=Q*C(:,:,k)*Q';B=B-diag(diag(B));D=D+sum(B(:).^2);
 end
 J=K-rho*D;
end
function Q=givensQ(t)
 a=t(1);b=t(2);c=t(3);
 G23a=[1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
 G12b=[cos(b) -sin(b) 0;sin(b) cos(b) 0;0 0 1];
 G23c=[1 0 0;0 cos(c) -sin(c);0 sin(c) cos(c)];
 Q=G23c*G12b*G23a;
end
function [meanCorr,amari]=metrics(y,s,G)
 cc=abs(corrcoef([y' s']));cc=cc(1:3,4:6);
 pp=perms(1:3);score=zeros(size(pp,1),1);
 for p=1:size(pp,1)
  score(p)=mean(cc(sub2ind([3 3],1:3,pp(p,:))));
 end
 meanCorr=max(score);
 G=abs(G);G=max(G,eps);
 rowSum=sum(bsxfun(@rdivide,G,max(G,[],2)),2)-1;
 colSum=sum(bsxfun(@rdivide,G,max(G,[],1)),1)-1;
 amari=(sum(rowSum)+sum(colSum))/(2*3*(3-1));
end
