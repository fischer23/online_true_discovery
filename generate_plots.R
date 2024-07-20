#This file creates the plots of the paper "Online closed testing with e-values"

#Load the simulated data that was generated with the file "generate_data.R"

load("results/TD_data.rda")
load("results/TD_data_taus.rda")

library(ggplot2)

col=c( "aquamarine", "purple", "cornflowerblue",  "pink", "orange", "red", "limegreen", "grey")


#Creates Figure 1 of the paper "Online true discovery guarantee with e-values"

plot_taus=ggplot(taus_df, aes(x=index, y=taus)) + 
  geom_line() +
  geom_line(aes(y=pi_A), col="grey", linetype = "dashed") +
  facet_grid(mu_A ~ pi_A_name)+
  xlab("Hypothesis index")+
  ylab("Estimated proportion of false hypotheses")+
  scale_x_continuous(breaks = seq(10, 90, 20), limits=c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) 

ggsave("results/plot_taus.pdf", plot=plot_taus, width=12, height=8)


#Creates Figure 2 of the paper "Online true discovery guarantee with e-values"

df_KR=TDP_df[TDP_df$Procedure %in% c("online-simple", "closed online-simple", "admissible online-simple", "true proportion"), ]
col_KR=c("aquamarine", "purple", "cornflowerblue", "grey")

plot_KR=ggplot(df_KR, aes(x=index, y=TDP)) + 
  geom_line(aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_KR)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) 

ggsave("results/plot_KR.pdf", plot=plot_KR, width=12, height=8)


#Creates Figure 3 of the paper "Online true discovery guarantee with e-values"

df_GRO=TDP_df[TDP_df$Procedure %in% c("GRO", "hedged GRO", "boosted GRO", "true proportion"), ]
col_GRO=c("pink", "orange", "red", "grey")

plot_GRO=ggplot(df_GRO, aes(x=index, y=TDP)) + 
  geom_line(aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_GRO)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("results/plot_GRO.pdf", plot=plot_GRO, width=12, height=8)


#Creates Figure 4 of the paper "Online true discovery guarantee with e-values"

df_best=TDP_df[TDP_df$Procedure %in% c("admissible online-simple", "boosted GRO", "calibrated", "true proportion"), ]
col_best=c("cornflowerblue", "red", "limegreen", "grey")

plot_best=ggplot(df_best, aes(x=index, y=TDP)) + 
  geom_line(aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_best)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("results/plot_best.pdf", plot=plot_best, width=12, height=8)