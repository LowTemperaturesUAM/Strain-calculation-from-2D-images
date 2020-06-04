%% 1. Load topography, plot topo & FFT
% MatrizTopo should be a square matrix

MatrizTopo=[1 0;0 1]; % <<<<< This is an input

figure(1)
imagesc(MatrizTopo)
axis square

FMatrizTopoOriginal=fft2(MatrizTopo);
FMatrizTopo=abs(fftshift(FMatrizTopoOriginal));
figure(2)
imagesc(FMatrizTopo)
axis square


%% 2. Calculate Q vectors from FFT
% Go to the FFT in figure 2 and locate the Bragg peaks with the data
% pointer. Also estimate the spatial extension of the peak (Lambda)
tam=size(MatrizTopo);
pix=tam(1);
qx=[1 0]; % <<<<< This is an input, in pixel units
qy=[0 1]; % <<<<< This is an input, in pixel units
Zero=[floor(pix/2)+1,floor(pix/2)+1];


Qx=(2*pi/pix)*(qx-Zero); % In reciprocal units
Qy=(2*pi/pix)*(qy-Zero); % In reciprocal units
ax=(2*pi/norm(Qx)); % Lattice parameter in x, in pixels
ay=(2*pi/norm(Qy)); % Lattice parameter in y, in pixels

%% 3. Calculate Tx(x,y) and  Ty(x,y)
% They will be complex matrices containing the information about the phase
% lag in the topography with respect to the perfect lattice formed by the
% vectors Qx and Qy
im=sqrt(-1);
Lambda_pix=1; % <<<<< This is an input, in pixel units
Lambda=(2*pi/pix)*Lambda_pix; % in reciprocal space, in reciprocal pixel units
OneOverLambda=2*pi/Lambda; % in real space, in pixels

Tx=zeros(pix);
Ty=zeros(pix);

% This loop might take several hours to run, depending on MatrizTopo size
% and the computing power
for i=1:pix
    for j=1:pix
        cont=1;
        r=[j,i];
        for k=1:pix
            for l=1:pix
                rprima=[l,k];
                sumX(cont)=MatrizTopo(k,l).*exp(-im*dot(Qx,rprima)).*((Lambda^2/(2*pi)).*exp(-0.5*Lambda^2*norm(r-rprima)^2));
                sumY(cont)=MatrizTopo(k,l).*exp(-im*dot(Qy,rprima)).*((Lambda^2/(2*pi)).*exp(-0.5*Lambda^2*norm(r-rprima)^2));
                cont=cont+1;
            end
        end
        Tx(i,j)=sum(sumX);
        Ty(i,j)=sum(sumY);
        clear sumX
        clear sumY
    end
end
%% 4. Calculate the angle lag matrices ATx, ATy
ATx=atan2(imag(Tx),real(Tx));
ATy=atan2(imag(Ty),real(Ty));

figure(21)
imagesc(ATx)
axis square

figure(22)
imagesc(ATy);
axis square

% We can visualize ATx and ATy and observe if there are any 2pi phase
% slips that need to be corrected

%% 5.1. Correct phase slips in ATx (if necessary)

Cx=ATx; % Dummy matrix to not modify the original ATx
delay=0; % number of pixels that the algorith skips before it starts 
         % correcting the phase slips. This might be necessary if there
         % any edge distortions in ATx.
tol=6; % any phase jump of 2pi will overcome this threshold
dCx=diff(Cx,1,1); % we locate the phase jumps by differentiating
MCx=zeros(pix);
for i=2:pix-1
    for j=2:pix-1
        MCx(i,j)=MCx(i-1,j);
        jump=dCx(i,j);
        if jump>tol
            MCx(i,j)=MCx(i-1,j)-1;
        elseif jump<-tol
            MCx(i,j)=MCx(i-1,j)+1;
        end
    end
end
MCx(2:pix-2*delay,:)=MCx(1:pix-2*delay-1,:);
% MC should be an "offset shadow" of ATx
figure(1401)
imagesc(MCx)
axis square

% We cut MC leaving the "delay" boder out
MCrecx=2*pi*MCx(delay:pix-delay-1,delay:pix-delay-1);
figure(1402)
imagesc(MCrecx)
axis square

Crecx=ATx(delay:pix-delay-1,delay:pix-delay-1);
Creccorrx=Crecx+MCrecx;

figure(1403)
imagesc(Creccorrx)
axis square

ATxCorr=Creccorrx;
Zx=zeros(pix);
Zx(delay:pix-delay-1,delay:pix-delay-1)=ATxCorr;
ATx=Zx; % Now we finally correct ATx

%% 5.2. Correct phase slips in ATy (if necessary)

Cy=ATy;
delay=0;
tol=6;
dCy=diff(Cy,1,1);
MCy=zeros(pix);
for i=2:pix-1
    for j=2:pix-1
        MCy(i,j)=MCy(i-1,j);
        jump=dCy(i,j);
        if jump>tol
            MCy(i,j)=MCy(i-1,j)-1;
        elseif jump<-tol
            MCy(i,j)=MCy(i-1,j)+1;
        end
    end
end
MCy(2:pix-2*delay,:)=MCy(1:pix-2*delay-1,:);
figure(1401)
imagesc(MCy)
axis square

MCrecy=2*pi*MCy(delay:pix-delay-1,delay:pix-delay-1);
figure(1402)
imagesc(MCrecy)
axis square

Crecy=ATy(delay:pix-delay-1,delay:pix-delay-1);
Creccorry=Crecy+MCrecy;

figure(1403)
imagesc(Creccorry)
axis square
ATyCorr=Creccorry;
Zy=zeros(pix);
Zy(delay:pix-delay-1,delay:pix-delay-1)=ATyCorr;
ATy=Zy;


%% 6. Calculate u(r)
% ux(x,y) and uy(x,y) can be calculated by solving a linear set of
% algebraic equations

matQ=[Qy;Qx];
sol=inv(matQ);
for i=1:pix
    for j=1:pix
        matAng=[-ATx(i,j),-ATy(i,j)];
        uy(i,j)=sol(1,1).*matAng(1)+sol(1,2)*matAng(2);
        ux(i,j)=sol(2,1).*matAng(1)+sol(2,2)*matAng(2);
    end
end

%% 7. Remove artificial plane from ux, uy
% Since Qx and Qy don't have infinite resolution in reciprocal space, the
% phase lag from the lattice might have a small offset phase that
% accumulates from pixel to pixel. This offset might also come from a
% uniform piezo or thermal drift.  We cancel this artifact linear term by
% substracting the global plane to ux and uy

borde=ceil(OneOverLambda/2); % distance from edge where the result is not 
                             % reliable due to poor statistics   
uxred=ux(borde+1:pix-borde,borde+1:pix-borde); % Crop the relevant central part of ux
uyred=uy(borde+1:pix-borde,borde+1:pix-borde); % Same in uy
tamred=pix-2*borde; % Size of the new reduced uxred & uyred

MeanXux=mean(uxred,1);
MeanYux=mean(uxred,2);
MeanXuy=mean(uyred,1);
MeanYuy=mean(uyred,2);

TiltXux=diff(MeanXux);
TiltYux=diff(MeanYux);
TiltXuy=diff(MeanXuy);
TiltYuy=diff(MeanYuy);

MeanTiltXux=mean(TiltXux);
MeanTiltYux=mean(TiltYux);
MeanTiltXuy=mean(TiltXuy);
MeanTiltYuy=mean(TiltYuy);

planeXux=MeanTiltXux.*(1:tamred);
planeYux=MeanTiltYux.*(1:tamred);
planeYux=planeYux';
planeXux=repmat(planeXux,tamred,1);
planeYux=repmat(planeYux,1,tamred);
planeux=planeXux+planeYux;
uxredplane=uxred-planeux;

planeXuy=MeanTiltXuy.*(1:tamred);
planeYuy=MeanTiltYuy.*(1:tamred);
planeYuy=planeYuy';
planeXuy=repmat(planeXuy,tamred,1);
planeYuy=repmat(planeYuy,1,tamred);
planeuy=planeXuy+planeYuy;
uyredplane=uyred-planeuy;

% Now we remove the offset in uxredplane, uyredplane
uxf=uxredplane-mean(mean(uxredplane));
uyf=uyredplane-mean(mean(uyredplane));
% We generate the vector field in the same window of size pix of
% MatrixzTopo
uxfc=zeros(pix);
uyfc=zeros(pix);
uxfc(borde+1:pix-borde,borde+1:pix-borde)=uxf;
uyfc(borde+1:pix-borde,borde+1:pix-borde)=uyf;

% We finally find our displacement field
% u(r) = [ uxfc(i,j) , uyfc(i,j) ]

%% 8. We calculate uxx, uyy and s
% By differentiating ux and uy we can access to the 2D strain tensor S
% components uxx and uyy.

% S(x,y) = [uxx(x,y) uxy(x,y) ; uyx(x,y) uyy(x,y)]

% The classic lattice strain s is the trace of S. The minus sign is
% arbitrary and it's sele
% s(x,y)=-0.5*(uxx+uyy)

duxf=diff(uxf,1,2);
uxx=-duxf(1:pix-2*borde-1,:);
duyf=diff(uyf,1,1);
uyy=-duyf(:,1:pix-2*borde-1);

s=zeros(pix);
s(borde+1:pix-borde-1,borde+1:pix-borde-1)=-0.5*(uxx+uyy);
% t(borde+1:pix-borde-1,borde+1:pix-borde-1)=-0.5*(uxx-uyy);

figure(10)
imagesc(s)
axis square

figure(4)
imagesc(uxfc)
axis square

figure(5)
imagesc(uyfc)
axis square

% figure(6)
% imagesc(uxx)
% axis square
% 
% figure(7)
% imagesc(uyy)
% axis square

% figure(5)
% imagesc(MatrizTopo), hold on
% colormap gray
% quiver(uxfc,uyfc,'r','linewidth',2)
% axis square


%% 9. We store everything in an InfoStruct for latter use
InfoStruct_Image.Lambda_pix=Lambda_pix;
InfoStruct_Image.Lambda=Lambda;
InfoStruct_Image.MatrizTopo=MatrizTopo;
InfoStruct_Image.FMatrizTopo=FMatrizTopo;
InfoStruct_Image.qx=qx;
InfoStruct_Image.qy=qy;
InfoStruct_Image.Cero=Zero;
InfoStruct_Image.Qx=Qx;
InfoStruct_Image.Qy=Qy;
InfoStruct_Image.ax=ax;
InfoStruct_Image.ay=ay;
InfoStruct_Image.ATx=ATx;
InfoStruct_Image.ATy=ATy;
InfoStruct_Image.MC=MC;
InfoStruct_Image.ATxCorr=ATxCorr;
InfoStruct_Image.ATyCorr=ATyCorr;
InfoStruct_Image.uxfc=uxfc;
InfoStruct_Image.uyfc=uyfc;
InfoStruct_Image.uxx=uxx;
InfoStruct_Image.uyy=uyy;
InfoStruct_Image.s=s;
