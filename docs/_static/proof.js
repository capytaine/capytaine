  $(document).ready(function() {
      $(".proof-type-proof > *").hide();
      $(".proof-type-proof .proof-title").show();
      $(".proof-type-proof .proof-title").click(function() {
          $(this).parent().children().not(".proof-title").toggle(400);
          $(this).parent().children(".proof-title").toggleClass("open");
      })
  });
