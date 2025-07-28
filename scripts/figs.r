library(dplyr)
library(ggplot2)
library(here)

load(file=here("results/run_sim.RData"))

l1 %>% 
    filter(out == "prs", family_rank == 1) %>%
    group_by(gen, what, h2, rho) %>%
    summarise(r=median(r^2)) %>%
    mutate(gen=gen+1, what = case_when(what == "rare" ~ "Null rare variants", what == "common" ~ "Null common variants", what == "prs" ~ "Causal variants")) %>%
    mutate(rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2lab ~ rholab, labeller = label_parsed) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and PRS"), colour="Variant") +
    scale_colour_brewer(type="qual") +
    theme_minimal()
ggsave(here("results/run_sim_prs.pdf"), width=7, height=7)

l1 %>% 
    filter(out == "y", family_rank == 1) %>%
    group_by(gen, what, h2, rho) %>%
    summarise(r=median(r^2)) %>%
    mutate(gen=gen+1, what = case_when(what == "rare" ~ "Null rare variants", what == "common" ~ "Null common variants", what == "prs" ~ "Causal variants")) %>%
    mutate(rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2lab ~ rholab, labeller = label_parsed) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and trait"), colour="Variant") +
    scale_colour_brewer(type="qual") +
    theme_minimal()
ggsave(here("results/run_sim_trait.pdf"), width=7, height=7)


load(file=here("results/run_sim2.RData"))
temp <- l2 %>% 
    filter(out == "y") %>%
    group_by(family_rank, what, h2, rho, gen) %>%
    summarise(r=median(r^2)) %>%
    mutate(
        rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)
    )

str(temp)
temp <- mutate(temp, 
    what = case_when(
        what == "rare" ~ "Null rare variants", 
        what == "common" ~ "Null common variants", 
        what == "prs" ~ "Causal variants"
    )
)

p1 <- temp %>%
    ggplot(., aes(x=family_rank, y=r)) + 
    geom_point(aes(colour=as.factor(gen))) +
    geom_smooth(aes(colour=as.factor(gen)), se=FALSE) +
    facet_grid(. ~ what) +
    labs(x="Family phenotype rank", y=expression("Median" ~ R^2 ~ "between variant and trait"), linetype="Variant", colour="Generation") +
    scale_colour_brewer(type="qual") +
    theme_minimal() +
    theme(axis.text.x= element_text(angle=90, hjust=0.5))
ggsave(p1, file="l2.pdf", width=14)
