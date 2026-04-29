# Source - https://stackoverflow.com/a/15473259
# Posted by user974514, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-29, License - CC BY-SA 3.0

dat1<-cbind(dat,aux=rep(1,length(dat[,1]))) 
dat1<-within(dat1, {aux = unlist(by(aux,Treatment,cumsum))})
dat1$aux<-dat1$aux+as.numeric(dat1$Treatment)*10

dat1

ggplot(dat1, aes(x=aux, y = Meas, ymin = Meas - SD/2, ymax = Meas + SD/2)) +
  geom_linerange(aes(color = Temp), size = 1, alpha = 0.5) +geom_point(aes(color = Temp, shape = Temp))+
  scale_x_continuous("Treatment",breaks=c(13.5,23.5), labels=c("A","B")) + # here you define coordinates for A and B 
  theme_bw()


################################################################################

dat1 <- read.table('/Users/songweizhi/Desktop/Sponge_r226/09_Dating/VisHPD95_20260415/dating_results_MCMCTree_vs_RelTime_table.txt', header = T)


label_str = 'LCA_AOA,LCA_f__Nitrosocaldaceae,LCA_f__UBA213,LCA_f__Nitrososphaeraceae,LCA_f__Nitrosopumilaceae,LCA_D1,LCA_D3,LCA_D4,LCA_D5,LCA_D6,LCA_D7,LCA_D8,LCA_D9,LCA_D10,LCA_D11'
label_list = unlist(strsplit(label_str, split = ","))

ggplot(dat1, aes(x=Post_X, y = Mean, ymin = Low, ymax = High)) +
  geom_linerange(aes(color = Color), size = 1, alpha = 0.5) +geom_point(aes(color = Color, shape = Shape))+
  scale_x_continuous("Treatment",breaks=c(14.5,69.5,118.5,176.5,234.5,292.5,350.5,408.5,466.5,524.5,582.5,640.5,698.5,756.5,814.5), labels=label_list) + # here you define coordinates for A and B 
  theme_bw()

