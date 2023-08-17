library(caret)
set.seed(2022-12-09)


control = trainControl(method="repeatedcv", 
                       number = 10, 
                       repeats = 10,
                       classProbs = TRUE,
                       savePredictions=TRUE)



mod_rf = train((xgTissue) ~ .,
               data = train.data_model, 
               method='rf', 
               trControl = control)

mod_svm = train(xgTissue ~ .,
               data = train.data_model, 
               method='svmLinear', 
               trControl = control)

mod_xgbTree = train(xgTissue ~ .,
                data = train.data_model, 
                method='xgbTree', verbosity = 0,
                trControl = control)

# mod_xgbLinear = train(xgTissue ~ .,
#                 data = train.data_model,
#                 method='xgbLinear',
#                 trControl = control)

explainer_rf  <- explain(mod_rf, label = "RF",
                         data = test.data_model[,2:ncol(test.data_model)], 
                         y = test.data_model$xgTissue %>% as.factor())

explainer_svm  <- explain(mod_svm, label = "SVM",
                          data = test.data_model[,2:ncol(test.data_model)], 
                          y = test.data_model$xgTissue %>% as.factor())

explainer_xgbtree  <- explain(mod_xgbTree, label = "xgbTree",
                          data = test.data_model[,2:ncol(test.data_model)], 
                          y = test.data_model$xgTissue %>% as.factor())

# explainer_xgbLinear  <- explain(mod_xgbLinear, label = "xgbLinear",
#                           data = test.data_model[,2:ncol(test.data_model)], 
#                           y = test.data_model$xgTissue %>% as.factor())

# explainer_knn  <- explain(mod_knn, label = "knn",
#                           data = test.data_model[,2:ncol(test.data_model)], 
#                           y = test.data_model$xgTissue %>% as.factor())
# 


pred <- predict(mod_xgbTree, train.data_model[,-1] %>% as.matrix(),  type = "prob") 
pred$max <- apply(pred, 1, max, na.rm=TRUE)
pred <- pred %>%  mutate('class'=names(.)[apply(., 1, which.max)]) 
test_pred <- bind_cols(test.data_model, pred)
# Plot values for prediction by tissue
ggplot(test_pred, aes(x=xgTissue, y=class)) + geom_jitter()



mp_rf  <- model_performance(explainer_rf)
mp_svm  <- model_performance(explainer_svm)
mp_xgbtree  <- model_performance(explainer_xgbtree)

mpart_rf <- model_parts(explainer_rf, type = 'raw')
mpart_svm <- model_parts(explainer_svm, type = 'raw')
mpart_xgbtree <- model_parts(explainer_xgbtree, type = 'raw')


plot(mpart_svm) + plot(mpart_rf) + plot(mpart_xgbtree)


svm <- svm(as.factor(xgTissue) ~ ., data=train.data_model, type = 'C-classification')