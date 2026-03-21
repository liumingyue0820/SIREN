# 1. Loading R packages ----
library(vip)
library(magrittr)
library(tidyverse)
library(fastshap)
library(shapviz)

# 2. Evaluation of feature importance using VIP scores ----
vip_result <-vip(final_model_enrich)
vip_result$data

# 3. Evaluation of feature importance using SHAP values ----
final_matrix_enrich<-as.data.frame(final_matrix_enrich)
final_matrix_enrich[]<-lapply(final_matrix_enrich,as.numeric)
final_matrix_enrich<-as.matrix(final_matrix_enrich)
newdata<-final_matrix_enrich
matadata$irAE_grade_binary
model<-final_model_enrich

set.seed(123)
shap <- explain(model,X=newdata,nsim=10,
                pred_wrapper = function(model,newdata){
                  predict(model, newx = newdata, s = "lambda.1se", type = "response")[,1]
                }
)
shap_handle <- shap %>% as.data.frame() %>% mutate(id=1:n()) %>% pivot_longer(cols = -id,names_to = "feature", values_to = "shap") 
data2 <- newdata %>% as.data.frame() %>% mutate(id=1:n()) %>% pivot_longer(cols = -id,names_to = "feature",values_to = "value")
combined_data <- left_join(
  shap_handle,
  data2,
  by = c("id", "feature")  
)
feature_importance <- combined_data %>%
  group_by(feature) %>%                          
  summarise(
    mean_abs_shap = mean(abs(shap)),             
    mean_shap = mean(shap),                      
    freq_non_zero = sum(shap != 0) / n()         
  ) %>%
  arrange(desc(mean_abs_shap))
combined_data <- combined_data %>%
  mutate(
    feature = factor(
      feature,
      levels = feature_importance$feature  
    )
  )

# 4. Evaluation of collinearity using the Phi coefficient ----
phi <- function(x, y) {
  contingency_table <- table(x, y)
  (contingency_table[1,1] * contingency_table[2,2] - contingency_table[1,2] * contingency_table[2,1]) /
    sqrt(sum(contingency_table[,1]) * sum(contingency_table[,2]) * 
           sum(contingency_table[1,]) * sum(contingency_table[2,]))
}
gene_pairs_data<-as.data.frame(TrainData)
gene_pairs_data[]<-lapply(gene_pairs_data,as.numeric)

n_pairs <- ncol(gene_pairs_data)
phi_matrix <- matrix(NA, nrow = n_pairs, ncol = n_pairs)
p_value_matrix <- matrix(NA, nrow = n_pairs, ncol = n_pairs)
colnames(phi_matrix) <- colnames(gene_pairs_data)
rownames(phi_matrix) <- colnames(gene_pairs_data)
colnames(p_value_matrix) <- colnames(gene_pairs_data)
rownames(p_value_matrix) <- colnames(gene_pairs_data)

for (i in 1:n_pairs) {
  for (j in 1:n_pairs) {
    if (i != j) { 
      x <- gene_pairs_data[, i]
      y <- gene_pairs_data[, j]
      contingency_table <- table(x, y)
      phi_value <- phi(x, y)
      phi_matrix[i, j] <- phi_value
      chisq_result <- chisq.test(contingency_table, correct = FALSE)
      p_value_matrix[i, j] <- chisq_result$p.value
    }
  }
}
phi_matrix
p_value_matrix


