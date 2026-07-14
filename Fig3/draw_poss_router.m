
pathsa=('E:\资料\QRN\202407\20240725-wit\');
Fontsz=10;
Fontna='Times New Roman';
angl=0;
erwid=2;   %errorbar粗细
YName='Correlation'

load('E:\资料\QRN\202407\20240725-wit\CC_router.mat')
A(1,:)='A1';
A(2,:)='A2';
B(1,:)='B1';
B(2,:)='B2';
B(3,:)='B3';
C(1,:)='C1';
C(2,:)='C2';
C(3,:)='C3';
C(4,:)='C4';
D(1,:)='D1';
D(2,:)='D2';
D(3,:)='D3';
D(4,:)='D4';
n=1;
Yar=[B;C;D];

for ia=1:2                   %22
    for iaa=1:11
    Xsi(n,:)=[A(ia,:),Yar(iaa,:)];
    n=n+1;
    end
end
Yar=[C;D];
for ib=1:3                   %24 
    for ibb=1:8
    Xsi(n,:)=[B(ib,:),Yar(ibb,:)];
    n=n+1;
    end
end
Yar=D;
for ic=1:4                   %16
    for icc=1:4
        Xsi(n,:)=[C(ic,:),Yar(icc,:)];
        n=n+1;
    end
end

Ysi=zeros(62,1);
CYsi=zeros(62,1);

for j=1:8
    for i=1:62
        for z=1:28
            if all(Xsi(i,:)==YM(z,:,j)) || all(Xsi(i,:)== [YM(z,3:4,j),YM(z,1:2,j)])
                    Ysi(i)=dcf(j,z)+Ysi(i);
                    CYsi(i)=CYsi(i)+1;
            end
        end
    end
end

YYsi=Ysi./CYsi;
all_YYsi=sum(YYsi);
de_YYsi=YYsi.^0.5;
de_all=sum(YYsi)^0.5;
u=YYsi;
v=all_YYsi;
sigma=(de_YYsi.^2/(all_YYsi^2)+de_all^2*YYsi.^2/(all_YYsi^4)).^0.5;




Z_Y=YYsi/sum(YYsi);

xx=[1:62];
hold on
h1=bar(Z_Y);
h1_error=errorbar(xx,Z_Y,sigma);
theoypos=[6 11 12 34 45 47 52];
valthe=zeros(62,1);
valthe(theoypos)=[1/6 1/6 1/6 1/6 1/6 1/12 1/12];

h2=bar(-valthe);
yticks=[1 0.75 0.50 0.25 0 0.25 0.50 0.75 1];
ylim([-0.6,0.6])
Xsii=Xsi;
Xsii([2:61],:)=" ";
set(gca,'ytick',[-1:0.25:1],'yticklabels',yticks,'Fontname',Fontna,'Fontsize',Fontsz);
set(gca,'xtick',[1:62],'xticklabels',Xsi,'Fontname',Fontna,'Fontsize',Fontsz,'XTickLabelRotation',angl);
ylabel(YName);
%legend({'experiment','theory'});
set(gcf,"Position",[0 0 6000 300]);
%saveas(gca,[pathsa,['CC_All','Z']],'jpg');
Posb(1)=sum(YYsi([1:22]));
err_Pos(1)=sum(sigma([1:22]).^2);
Posb(2)=sum(YYsi([1:3 12:14 23:46]));
err_Pos(2)=sum(sigma([1:3 12:14 23:46]).^2);
Posb(3)=sum(YYsi([4:7 15:18 23:26 31:34 39:42 47:62]));
err_Pos(3)=sum(sigma([4:7 15:18 23:26 31:34 39:42 47:62]).^2);
Posb(4)=sum(YYsi([8:11 19:22 27:30 35:38 43:46 47:62]));
err_Pos(4)=sum(sigma([8:11 19:22 27:30 35:38 43:46 47:62]).^2);
figure
h4=bar(Posb/sum(Posb));
hold on
h4_err=errorbar([1:4],Posb/sum(Posb),err_Pos)
h5=plot([0,6],[0.25,0.25])
set(gca,'ytick',[-1:0.25:1],'yticklabels',yticks,'Fontname',Fontna,'Fontsize',Fontsz);
set(gca,'xtick',[1:4],'xticklabels',['A';'B';'C';'D'],'Fontname',Fontna,'Fontsize',Fontsz,'XTickLabelRotation',angl);
ylabel(YName);


