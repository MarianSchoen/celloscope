function getCookies(){
  var res = Cookies.get();
  Shiny.setInputValue('cookies', res);
}


$(document).on('shiny:connected', function(ev){
  // set a random number to identify this session: 
  Cookies.set("session-number", Math.floor(Math.random()* 100000)); 
  getCookies();
})