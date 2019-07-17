# COTSModel_Figure Creation

b = -0.4/22
yint = 1-b*3

df = data.frame(CCRatio = 1:35, Mortality = c(rep(1,2), yint+b*(3:25), rep(0.6, 10)))
                
p1 =ggplot(df, aes(x=CCRatio, y=Mortality)) + geom_line(color="coral", size=1.5) +
  theme_classic() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,35), expand = c(0,0)) +
  geom_vline(xintercept = c(3,25), linetype = "dashed") +
  labs(title="Ratio-Dependent COTS Mortality",
            x ="% Coral Cover / COTS (per Manta Tow)", y = "COTS Mortality") +
  annotate(geom="text", x=30, y=0.64, label="M = 0.6")

# Fecundity             
b = 0.8/25
yint = 0.2

df = data.frame(CCRatio = 0:35, Mortality = c(yint+b*(0:25), rep(1, 10)))

p2 = ggplot(df, aes(x=CCRatio, y=Mortality)) + geom_line(color="forestgreen", size=1.5) +
  theme_classic() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,35), expand = c(0,0)) +
  geom_vline(xintercept = c(25), linetype = "dashed") +
  labs(title="Ratio-Dependent COTS Fecundity",
       x ="% Coral Cover / COTS (per Manta Tow)", y = "COTS Fecundity") +
  annotate(geom="text", x=30, y=0.94, label="MaxFecund = 2e^7")

gridExtra::grid.arrange(p1,p2, ncol=2)
