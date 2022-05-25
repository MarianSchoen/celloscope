$(document).ready(function(){
  $('#dragdropupload_progress').on("DOMSubtreeModified",function(){
    
    var target = $('#dragdropupload_progress').children()[0];
    if(target.innerHTML === "Upload complete"){
      console.log('Change')
      target.innerHTML = '';      
    }
  });
});